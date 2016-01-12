#Crea un master frame usando como referencia(s) epocas dadas por el usuario
import os
import numpy as np
import glob
import sys
from astropy import wcs
from astropy.utils.console import ProgressBar, color_print
from astropy.io import fits
from scipy.optimize import curve_fit
from subprocess import check_output
from MP2015utils import *

#Para ignorar los warnings de WCS
import warnings
warnings.filterwarnings("ignore")

# PARAMETROS

iteraciones = 4     #Numero de iteraciones
threshold   = .5    #Porcentaje (0-1) en que debe aparecer una estrella para ser considerada
match_tol   = .3    #Tolerancia en arcsec para hacer los matches (0.3 arcsec ~ 1 pix en VVV)

# CARPETAS
radec_folder  = 'RADEC'  #Carpeta donde guarda los archivos con RA DEC
match_folder  = 'MATCH_MF'  #Carpeta donde guarda los matches
stilts_folder = os.path.dirname(os.path.realpath(__file__)) #STILTS debe estar con el archivo .py


#####Informacion extra
if len(sys.argv) == 1:

    print
    ee = color_print('Como ejecutar','yellow')
    print 'python', sys.argv[0], 'path/to/catalogs/'
    print
    print 'Archivo "\033[1;33mzinfo_img\033[0m" debe estar en la carpeta de los datos'
    print 'Se requiere un archivo "\033[1;33mref_master\033[0m" que contenga la(s) imagen(es) de referencia'
    print 'Outputs:'

    sys.exit(1)

#Ruta de los .dao, .fits y zinfo
data_folder = sys.argv[1]
if not data_folder.endswith('/'):
    data_folder += '/'

#Elimina los archivos que utiliza, para no confundirse con resultados anteriores
os.system('rm -r RADEC')
os.system('rm -r MATCH_MF')

#Crea las carpetas para guardar RADEC y matches
makedir(radec_folder)
makedir(match_folder)

#Lee ref_master y zinfo
color_print('\nLeyendo archivo con referencias...', 'cyan')
epn         = np.genfromtxt(data_folder+'zinfo_img', usecols=(0,), dtype='string')
se, el, yr  = np.genfromtxt(data_folder+'zinfo_img', unpack=True, usecols=(4,5,6))

km = np.array(['k' in e for e in epn])

refs = np.genfromtxt('ref_master', dtype='string')
kr   = refs[np.array(['k' in r for r in refs])]
print '\t%d referencias a usar' % len(kr)

#Busca archivos
color_print('\nBuscando archivos...', 'cyan')
k_files = np.sort(glob.glob('%s*k*.dao' % data_folder))
k_files = np.array([s for s in k_files if any(xs in s for xs in kr)])

print '\t%d de %d Ks encontradas' % (len(k_files), len(kr))

#Convierte a RADEC
color_print('\nConvirtiendo XY a RADEC...','cyan')
print '\tConvirtiendo Ks'
ProgressBar.map(XYtoRADEC,k_files,multiprocess=True)

'''
#Crea el CMD
color_print('Creando CMDs...', 'cyan')
jcmd = np.repeat(['RADEC/%s' % jr[0].replace('.dao', '.dat')], len(kr))
kcmd = np.array(['RADEC/%s' % k.replace('.dao', '.dat') for k in kr])
acmd = np.transpose([kcmd, jcmd])

ProgressBar.map(make_CMD, acmd, multiprocess=True)
'''

#Selecciona la epoca con mas estrellas como referencia
color_print('\nSeleccionando epoca de referencia...','cyan')
def wc(ep):
    out = check_output(['wc','-l',ep]).split()[0]
    return int(out)

nro_estrellas = np.zeros(len(kr))
with ProgressBar(len(kr)) as bar:
    for i in range(len(kr)):
        nro_estrellas[i] = wc(k_files[i])

referencia = k_files[np.argmax(nro_estrellas)]
print '\tReferencia seleccionada: %s' % referencia

#Matches de referencia con epocas
color_print('\nIniciando matches de cada epoca con la referencia...','cyan')
def match_epo(ep, cmd, ofn):
    com_epo = 'java -jar %s/stilts.jar tmatch2 ifmt1=ascii ifmt2=ascii matcher=sky ofmt=ascii values1="RA DEC" values2="RA DEC" \
            in1=%s in2=%s out=%s params=%.1f progress=none join=1and2' % (stilts_folder,cmd,ep,ofn,match_tol)

    os.system(com_epo)

krd    = glob.glob('RADEC/*k*.dat')
datref = glob.glob('RADEC/%s' % referencia.split('/')[-1].replace('.dao','.dat'))[0]
suf    = referencia.split('_')[-1].replace('.dao','.match')
'''
def me2(ep):
    ofn = ep.split('.')[0].replace('RADEC', 'MATCH') + '_' + suf
    match_epo(ep, datref, ofn)

ProgressBar.map(me2, krd, multiprocess=True)
'''

com_epo  = 'java -jar %s/stilts.jar tmatchn multimode=group nin=%d matcher=sky params=%.3f ' % (stilts_folder, len(kr), match_tol)
for i in range(1,len(kr)+1):
    com_epo += 'in%d=%s ifmt%d=ascii values%d="RA DEC" join%d=always ' % (i, krd[i-1], i, i, i)
com_epo += 'ocmd="addcol -before \$1 ID index" out=match_all.dat ofmt=ascii'
os.system(com_epo)

#Transformaciones lineales
color_print('Creando Master Frame...','cyan')

dataa = np.genfromtxt('match_all.dat', unpack=True)
equis = dataa[4::6].T

ma = [np.isfinite(x) for x in equis]
ma = np.sum(ma, axis=1) > 1

equis = np.nanmean(dataa[4::6].T, axis=1)[ma]
yes   = np.nanmean(dataa[5::6].T, axis=1)[ma]
ra    = np.nanmean(dataa[2::6].T, axis=1)[ma]
dec   = np.nanmean(dataa[3::6].T, axis=1)[ma]
mag   = np.nanmean(dataa[6::6].T, axis=1)[ma]
ids   = np.arange(len(equis))

fmt = '%d %f %f %.3f %.3f %.3f'
np.savetxt('MasterFrame.dat', np.transpose([ids, ra, dec, equis, yes, mag]), fmt=fmt)
