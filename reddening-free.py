import os
import sys
import numpy as np
from astropy.io import fits as pf
from astropy import wcs
from astropy.utils.console import ProgressBar, color_print

import warnings
warnings.filterwarnings("ignore")


folder = sys.argv[1]
zinfo  = folder + 'zinfo_img'

stilts_folder = os.path.dirname(os.path.realpath(__file__))

def makedir(directory):
	if not os.path.exists(directory):
		os.makedirs(directory)

def XYtoRADEC(cat):
	ffn = cat
    cfn = cat.replace('fits','dao')

	ids,x,y,mag = np.loadtxt(cfn,usecols=range(4),skiprows=3,unpack=True)
	ids = ids.astype(int)

	hdr    = pf.getheader(ffn)
	w      = wcs.WCS(hdr)
	ra,dec = np.transpose(w.wcs_pix2world(np.transpose([x,y]),1))

	head = 'ID RA DEC X Y MAG'
	fmt  = '%d %.7f %.7f %.3f %.3f %.3f'
    ofn  = './CMD/' + cfn.replace('dao','dat').split('/')[-1]
	np.savetxt(ofn, np.transpose([ids,ra,dec,x,y,mag]),header=head,fmt=fmt)

makedir('CMD')

#Lee el zinfo
files  = np.genfromtxt(zinfo,unpack=True,usecols=(0,),dtype='string')
seeing = np.genfromtxt(zinfo,unpack=True,usecols=(4,))

k_ma = np.array(['k' in f for f in files])
j_ma = np.array(['j' in f for f in files])
z_ma = np.array(['z' in f for f in files])
h_ma = np.array(['h' in f for f in files])
y_ma = np.array(['j' in f for f in files])

sk,sj,sz,sh,sy = seeing[k_ma],seeing[j_ma],seeing[z_ma],seeing[h_ma],seeing[y_ma]
k_file = files[k_ma][np.argmin(sk)]
j_file = files[j_ma][np.argmin(sj)]
z_file = files[z_ma][np.argmin(sz)]
h_file = files[h_ma][np.argmin(sh)]
y_file = files[y_ma][np.argmin(sy)]

#RADEC
go_files = [k_file,j_file,z_file,h_file,y_file]
ProgressBar.map(XYtoRADEC,go_files,multiprocess=True)

#Match
ejec = 'java -jar %s/stilts.jar tmatchn multimode=pairs nin=5 matcher=sky ' % stilts_folder
for i,f in enumerate(go_files)
    f = f.replace('.fits','dat')
    ejec += 'in%d=%s ifmt%d=ascii values%d="RA DEC" ' % (i+1, f, i+1, i+1)
ejec += 'join1=match out=./CMD/match.dat ofmt=ascii'
os.system(ejec)

#Indices
data = np.genfromtxt('./CMD/match.dat',unpack=True)

m1 = K - 1.08 * (Z - Y)
m2 = H - 1.13 * (J - K)
m3 = J - 1.03 * (Y - K)
m4 = K - 1.22 * (J - H)

c1 = (Y - J) - 1.14 * (J - H)
c2 = (Z - Y) - 0.99 * (Y - J)
c3 = (J - H) - 1.47 * (H - K)
c4 = (J - K) - 1.50 * (Z - Y)
