import numpy,pyfits,sys,os,subprocess,math

from tractor import *
from tractor.sdss import *
from astrometry.sdss import *
from astrometry.util.sdss_radec_to_rcf import *

sdss = DR9(basedir='data/unzip')

window_flist = 'window_flist.fits'

# Your target object, in decimal degrees:
ra,dec = 15.*(9.+51./60. + 22.57/3600.),(26.+35./60.+13.9/3600.)

# ra,dec = 123.523,58.478

# List of SDSS fields (run,camcol,field,ra,dec):
rcf = radec_to_sdss_rcf(ra, dec, tablefn=window_flist)

# The SDSS bands to read:
bands = 'griz'

# The RA,Dec box (ie, ra,dec plus margin in degrees)
m = 10./3600.
rlo = ra - m
rhi = ra + m
dlo = dec - m
dhi = dec + m

# Save all the Tractor Image objects here:
tims = []
# Read each SDSS field 
for (r,c,f,rr,dd) in rcf:

     for band in bands:
         tim,inf = get_tractor_image_dr9(
                 r, c, f, band, sdss=sdss,
                 nanomaggies=True, zrange=[-2,5],
                 roiradecbox=[rlo,rhi,dlo,dhi],
                 invvarIgnoresSourceFlux=True)
         if tim is None:
             continue
         print 'Read', tim.shape, 'pixels for', (r,c,f,band)
         tims.append(tim)


