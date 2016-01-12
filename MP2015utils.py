import numpy as np
from astropy import wcs
from astropy.io import fits
import os

stilts_folder = os.path.dirname(os.path.realpath(__file__))

def make_CMD(epx, tol=0.3):
    epk, epj = epx
    kn = int(epk.split('.')[0].split('-')[-1])
    jn = int(epj.split('.')[0].split('-')[-1])

    output = 'CMD_%03d-%03d.dat' % (kn, jn)

    #print '\tCreando %s' % output
    os.system('java -jar %s/stilts.jar tmatch2 \
               in1=%s values1="RA DEC" ifmt1=ascii \
               in2=%s values2="RA DEC" ifmt2=ascii \
               matcher=sky params=%f find=best join=1and2 \
               out=%s ofmt=ascii progress=none' % (stilts_folder,epk,epj,tol,output) )

def XYtoRADEC(ep, folder='RADEC'):
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

    np.savetxt('./%s/%s' % (folder,of),np.transpose([rid,ra,dec,rx,ry,rmag]),header=head,fmt=fmt)

def makedir(directory):
	if not os.path.exists(directory):
		os.makedirs(directory)

#Transformacion lineal
def linear(coords,a,b,c):
	x,y = coords
	return a + b*x + c*y
