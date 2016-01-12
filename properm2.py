import os
import sys
import glob
import numpy as np
from sklearn.neighbors import NearestNeighbors as NN
from astropy.utils.console import ProgressBar, color_print
from scipy.optimize import curve_fit
from fit import polyfitr

#PARAMETROS
vecinos = 56
output  = 'PM_'
locales = 'brgb.dat'
master  = 'Master.dat'
ma1,ma2 = 9, 20
ml1,ml2 = 11, 14
thresh  = 40

vc_pix = 0.34

match_folder = 'matched_epochs'
match_master = 'matched_mf'
pm_folder    = 'epochs_pm'

stilts_folder = os.path.dirname(os.path.realpath(__file__)) #STILTS debe estar con el archivo .py

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

#PIPELINE
folder = sys.argv[1]
makedir(match_folder)
makedir(match_master)
makedir(pm_folder)

color_print('Leyendo informacion de epocas','cyan')
se,el,yr = np.genfromtxt(folder+'zinfo_img',unpack=True,usecols=(4,5,6))
name     = np.genfromtxt(folder+'zinfo_img',unpack=True,usecols=(0,),dtype='string')
k_mask   = np.array(['k' in f for f in name])
sek,elk,yrk = np.transpose([se,el,yr])[k_mask].T
yrk = (yrk-yrk[0])/365.242199

color_print('Recopilando archivos de epocas...','cyan')
epochs = glob.glob('./%s/*.*' % match_folder)

color_print('Realizando match de la MF con las epocas','cyan')
ejecuta = 'java -jar %s/stilts.jar tmatch2 in1=./%s values1="ID" ifmt1=ascii ' % (stilts_folder, master)

def mf_match(ep):
    ej2 = 'in2=%s values2="ID_1" ifmt2=ascii icmd2=\'keepcols "ID_1 X Y"\' matcher=exact find=best join=1and2 out=./%s/%s ofmt=ascii progress=none ocmd="delcols ID_1"' % (ep, match_master, ep.split('/')[-1].replace('.match','.mfma'))
    os.system(ejecuta + ej2)

ProgressBar.map(mf_match,epochs,multiprocess=True)

color_print('Realizando transformaciones lineales','cyan')
matches = glob.glob('./%s/*.*' % match_master)
bid     = np.genfromtxt(locales,unpack=True,usecols=(0,))

def shift(ep):
    ids,x1,y1,mk,mj,x2,y2 = np.genfromtxt(ep,unpack=True)

    local_mask = np.in1d(ids,bid)
    lid,lx1,ly1,lx2,ly2  = np.transpose([ids,x1,y1,x2,y2])[local_mask].T

    loc_xy = np.transpose([lx2,ly2])
    nbrs   = NN(n_neighbors=vecinos, algorithm='auto').fit(loc_xy)

    coo_xy = np.transpose([x2,y2])
    dist, idx = nbrs.kneighbors(coo_xy)
    idx  = idx[:,1:]
    dist = dist[:,1:]

    ctx = np.zeros(x1.size)
    cty = np.zeros(y1.size)

    for i in range(x1.size):
        star_locales = loc_xy[idx[i]]
        ep1_x = lx1[idx[i]]
        ep1_y = ly1[idx[i]]

        poptx, pcovx = curve_fit(linear,star_locales.T,ep1_x)
        popty, pcovy = curve_fit(linear,star_locales.T,ep1_y)

        ctx[i] += linear([x2[i],y2[i]],*poptx)
        cty[i] += linear([x2[i],y2[i]],*popty)

    shift_x = x1 - ctx
    shift_y = y1 - cty

    hdr  = 'ID X Y MAG_K MAG_J PMX PMY'
    fmt  = '%d %.3f %.3f %.3f %.3f %f %f'
    data = np.transpose([ids,x1,y1,mk,mj,shift_x,shift_y])
    np.savetxt('./%s/%s' % (pm_folder,ep.split('/')[-1].replace('.mfma','.pm')), data, header=hdr, fmt=fmt)

ProgressBar.map(shift,matches,multiprocess=True)

color_print('Calculando PMs...','cyan')
shifts = np.sort(glob.glob('./%s/*k*.*' % pm_folder))
print '\tRealizando match de los desplazamientos por epoca...'
ejec  = 'java -jar %s/stilts.jar tmatchn multimode=pairs nin=%d matcher=exact ' % (stilts_folder, len(shifts)+1)
ejec += 'in1=%s ifmt1=ascii values1="ID" ' % master

for i in range(1,len(shifts)+1):
	ejec += 'in%d=%s ifmt%d=ascii values%d="ID" ' % (i+1, shifts[i-1], i+1, i+1)
ejec += ' join1=match out=./%s/match.pm ofmt=ascii' % (pm_folder)
os.system(ejec)

shift_data = np.genfromtxt('./%s/match.pm' % pm_folder, unpack=True)

ids,x,y,mag,col = shift_data[0:5]
dx  = shift_data[10::7]
dy  = shift_data[11::7]

PM_X, PM_Y = np.empty(dx.shape[1]), np.empty(dy.shape[1])
PM_X[:], PM_Y[:] = np.nan, np.nan

dx_fin = np.isfinite(dx)
dy_fin = np.isfinite(dy)

with ProgressBar(len(dx.T)) as bar:
    for i in range(len(dx.T)):
        if dx_fin[:,i].sum() > thresh:
            ma = dx_fin[:,i]
            xx = yrk[ma]
            yy = dx[:,i][ma]
            try:
                #coeffx  = polyfitr(xx,yy,order=1,clip=3)[0]
                coeffx  = np.polyfit(xx,yy,1)
		PM_X[i] = coeffx[0]

                yy = dy[:,i][ma]
		coeffx  = np.polyfit(xx,yy,1)
                #coeffy  = polyfitr(xx,yy,order=1,clip=3)[0]
                PM_Y[i] = coeffy[0]
            except:
                continue
        bar.update()

fmt = '%d %.3f %.3f %.3f %.3f %f %f'
hdr = 'ID X Y MAG_K COL_JK PM_X PM_Y'
np.savetxt(output+locales, np.transpose([ids,x,y,mag,col,PM_X,PM_Y]), fmt=fmt, header=hdr)

#10 7
