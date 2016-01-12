import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy import wcs
from astropy.utils.console import ProgressBar, color_print
from astropy.io import fits
from scipy.optimize import curve_fit

import warnings
warnings.filterwarnings("ignore")

#PARAMETROS
iteraciones = 4     #Numero de iteraciones antes de promediar
threshold   = 50    #Numero de veces que debe aparecer una estrella para ser considerada
match_tol   = 0.3   #Tolerancia en arcosegundos para el match (0.3 ~ 1 pix)

radec_folder  = 'RADEC'		 #Carpeta donde se guardan los archivos con RA y DEC
match_folder  = 'matched_epochs' #Carpeta con los matches de epocas con la de referencia (sin / al final)
stilts_folder = os.path.dirname(os.path.realpath(__file__)) #STILTS debe estar con el archivo .py
cmd_out       = 'CMD.dat'   #Archivo con el CMD (K ref + J ref)

#FUNCIONES
def makedir(directory):
	if not os.path.exists(directory):
		os.makedirs(directory)

def linear(coords,a,b,c):
	x,y = coords
	return a + b*x + c*y

def quad(coords,a,b,c,d,e,f):
	x,y = coords
	return a + b*x + c*y + d*np.square(x) + e*np.multiply(x,y) + f*np.square(y)

def XYtoRADEC(ep):
    cfn = ep
    ffn = ep.replace('.dao','.fits')

    rid,rx,ry,rmag = np.genfromtxt(cfn,usecols=range(4),skip_header=3,unpack=True)

    hdr    = fits.getheader(ffn)
    w      = wcs.WCS(hdr)
    rcoo   = np.transpose([rx,ry])
    ra,dec = np.transpose(w.wcs_pix2world(rcoo,1))

    head = 'ID RA DEC X Y MAG'
    fmt  = '%d %f %f %.3f %.3f %.3f'
    of   = ep.split('/')[-1].replace('.dao','.dat')
    np.savetxt('./%s/%s' % (radec_folder,of),np.transpose([rid,ra,dec,rx,ry,rmag]),header=head,fmt=fmt)


#PIPELINE
data_folder = sys.argv[1]
if not data_folder.endswith('/'):
    data_folder += '/'

makedir(radec_folder)
makedir(match_folder)

#Busca archivos .dao en la carpeta
color_print('Buscando archivos .dao...','cyan')
k_files = np.sort(glob.glob('%s*k*.dao' % data_folder))
j_files = np.sort(glob.glob('%s*j*.dao' % data_folder))
print '\tEncontradas %d epocas en Ks' % len(k_files)
print '\tEncontradas %d epocas en J' % len(j_files)

#Carga la info de seeing, anios, etc
color_print('Obteniendo informacion de seeing y BJD...','cyan')
eps      = np.genfromtxt(data_folder+'zinfo_img',unpack=True,usecols=(0,),dtype='string')
se,el,yr = np.genfromtxt(data_folder+'zinfo_img',unpack=True,usecols=(4,5,6),dtype=np.float64)

info_k    = np.char.find(eps,'k') > -1
info_j    = np.char.find(eps,'j') > -1

keps = eps[info_k]
jeps = eps[info_j]

kse,kel,kyr = np.transpose([se,el,yr])[info_k].T
jse,jel,jyr = np.transpose([se,el,yr])[info_j].T

msk_idx = np.argmin(kse)
msj_idx = np.argmin(jse)

print '\tSeeing minimo en Ks encontrado en %s (%f)' % (keps[msk_idx],kse[msk_idx])
print '\tSeeing minimo en J encontrado en %s (%f)' % (jeps[msj_idx],jse[msj_idx])

#Obtiene los RA-DEC
color_print('Convirtiendo XY a RADEC...','cyan')
print '\tConvirtiendo Ks'
ProgressBar.map(XYtoRADEC,k_files,multiprocess=True)
print '\tConvirtiendo J'
ProgressBar.map(XYtoRADEC,j_files,multiprocess=True)

#Crea el CMD
color_print('Creando CMD...','cyan')
K_ref = './%s/%s' % (radec_folder,k_files[msk_idx].split('/')[-1].replace('.dao','.dat'))
J_ref = './%s/%s' % (radec_folder,j_files[msj_idx].split('/')[-1].replace('.dao','.dat'))

com_cmd = 'java -jar %s/stilts.jar tmatch2 ifmt1=ascii ifmt2=ascii matcher=sky ofmt=ascii values1="RA DEC" values2="RA DEC" \
        in1=%s in2=%s out=%s params=%.1f progress=none join=all1' % (stilts_folder,K_ref,J_ref,cmd_out,match_tol)
os.system(com_cmd)

#Match de epocas con la de referencia (CMD)
makedir(match_folder)

def match_epo(ep):
    com_epo = 'java -jar %s/stilts.jar tmatch2 ifmt1=ascii ifmt2=ascii matcher=sky ofmt=ascii values1="RA_1 DEC_1" values2="RA DEC" \
            in1=%s in2=%s out=./%s/%s params=%.1f progress=none join=1and2' % (stilts_folder,cmd_out,ep,match_folder,ep.split('/')[-1].replace('.dat','.match'),match_tol)

    os.system(com_epo)

K_dat = np.sort(glob.glob('./%s/*k*.dat' % radec_folder))
print '\tRealizando match de %d epocas Ks con la de referencia' % len(K_dat)
ProgressBar.map(match_epo,K_dat,multiprocess=True)

#MasterFrame
color_print('Creando Master Frame...','cyan')
print '\tNumero de iteraciones: %d' % iteraciones

mid,mx,my,mk,mj = np.genfromtxt('CMD.dat',usecols=(0,3,4,5,11),unpack=True)
K_matches = np.sort(glob.glob('./%s/*k*.match' % match_folder))

for i in range(iteraciones):
    counts = np.zeros((len(K_matches),mid.size))
    par_x  = np.zeros((len(K_matches),mid.size))
    par_y  = np.zeros((len(K_matches),mid.size))

    with ProgressBar(len(K_matches)) as bar:
	    for i,ep in enumerate(K_matches):
        	eid,ex,ey = np.genfromtxt(ep,unpack=True,usecols=(0,16,17))
        	mainep    = np.in1d(mid,eid)

        	ep_coords = np.transpose([ex,ey])
        	mmx, mmy  = mx[mainep], my[mainep]

        	poptx, pcovx = curve_fit(linear,ep_coords.T,mmx)
        	popty, pcovy = curve_fit(linear,ep_coords.T,mmy)

        	ctx, cty = linear(ep_coords.T,*poptx), linear(ep_coords.T,*popty)

        	counts[i][mainep] += 1
        	par_x[i][mainep]  += ctx
        	par_y[i][mainep]  += cty

		bar.update()

    id_counts = np.sum(counts,axis=0)   #Nro de veces que cada estrella fue encontrada
    ep_counts = np.sum(counts,axis=1)   #Nro de estrellas en cada epoca

    sum_x = np.sum(par_x,axis=0)
    sum_y = np.sum(par_y,axis=0)

    print id_counts
    print np.nanmax(par_x),np.nanmin(par_x), np.nanmean(par_x)
    print np.nanmax(par_y),np.nanmin(par_y), np.nanmean(par_y)

    mx,my = np.divide(sum_x,id_counts), np.divide(sum_y,id_counts)

color_print('Guardando MasterFrame','lightcyan')
thr_idx     = id_counts >= threshold
masterframe = np.transpose([mid,mx,my,mk,mj-mk])[thr_idx]

header = 'ID X Y MAG_K COL_JK'
fmt    = '%d %.3f %.3f %.3f %.3f'

np.savetxt('Master.dat',masterframe,fmt=fmt,header=header)

#GRAFICOS
color_print('Generando graficos...','cyan')
ep_stars = np.zeros(fits_k.size)

sys.exit()

import subprocess as sp
for i in range(K_dat.size):
	ep_stars[i] += (int(sp.check_output('wc -l %s'%(folder+fits_k[i].replace('fits','dao')), stderr=sp.STDOUT, shell=True).split()[0])-3)

import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=[6*2,4.5*3],nrows=3)

ax[0].fill_between(epocas,0,ep_stars,color='black',alpha=.66,zorder=-2)
ax[0].fill_between(epocas,0,ep_counts,color='orange',alpha=.66,zorder=0)
ax[0].plot(epocas,ep_stars,'-ok',alpha=.75,zorder=-1)
ax[0].plot(epocas,ep_counts,'-o',color='orange',alpha=.75,zorder=1)
ax[0].set_xlim(epocas.min()-1,epocas.max()+1)
ax[0].set_xlabel(r'Epocas')
ax[0].set_ylabel(r'Nro Estrellas')

ax[1].scatter(mx,id_counts,s=5,c=mk,edgecolor='',alpha=.5)
ax[1].set_xlabel(r'x')
ax[1].set_ylabel(r'Nro Epocas')
ax[1].set_xlim(-1,mx.max()+1)

ax[2].scatter(my,id_counts,s=5,c=mk,edgecolor='',alpha=.5)
ax[2].set_xlabel(r'y')
ax[2].set_ylabel(r'Nro Epocas')
ax[2].set_xlim(-1,mx.max()+1)

fig.savefig(results+'Master.png',dpi=300,bbox_inches='tight')
