'''
This file is part of the LensFinder project.
Copyright 2012 David W. Hogg (NYU) and Phil Marshall (Oxford).

Description
-----------

Running the tractor on PS1 *single object cutout* images.
Read in an image and its weight map, guess the PSF, put an object at the
centre image of the field and update the catalog.

Example use
-----------

 python ps1tractor.py -v \
    examples/H1413+117_10x10arcsec_55377.34051_z_sci.fits \
    examples/H1413+117_10x10arcsec_55377.34051_z_var.fits

'''

if __name__ == '__main__':
      import matplotlib
      matplotlib.use('Agg')

import os
import logging
import numpy as np
import pylab as plt
import pyfits

from astrometry.util.file import *
from astrometry.util.util import Tan
from astrometry.util.pyfits_utils import *

import tractor
from tractor import sdss_galaxy as gal
import lensfinder

# ============================================================================

def ps1tractor():

   from optparse import OptionParser
   import sys

   # Set available options:
   parser = OptionParser(usage=('%prog <sci> <var>'))
   # Verbosity:
   parser.add_option('-v', '--verbose', dest='verbose', action='count', \
                     default=False, help='Make more verbose')
   
   # Read in options and arguments - note only sci and wht images are supplied:
   opt,args = parser.parse_args()
   
   if len(args) != 2:
      parser.print_help()
      sys.exit(-1)
   scifile, varfile = args
 
   # -------------------------------------------------------------------------
   # Logging to terminal:
   
   if opt.verbose:
      lvl = logging.DEBUG
   else:
      lvl = logging.INFO
   logging.basicConfig(level=lvl, format='%(message)s', stream=sys.stdout)

   # -------------------------------------------------------------------------
   # Get sci and wht images, and make mask:
   
   hdulist = pyfits.open(scifile)
   sci = hdulist[0].data
   hdr = hdulist[0].header
   hdulist.close()
   NX,NY = sci.shape

   hdulist = pyfits.open(varfile)
   var = hdulist[0].data
   hdulist.close()
   assert sci.shape == var.shape
   
   mask = numpy.ones([NX,NY],dtype=np.int16)
   mask[numpy.where(var == 0)] = 0
   
   # Convert var to wht, and find median uncertainty as well:
   invvar = 1.0/var
   # Assign zero weight to var=nan, var<=0:
   invvar[var != var] = 0.0
   invvar[var <= 0] = 0.0
   
   good = np.where(invvar > 0)
   bad = np.where(invvar == 0)
   
   # Zero out sci image where wht is 0.0:
   sci[bad] = 0.0

   assert(all(np.isfinite(sci.ravel())))
   assert(all(np.isfinite(invvar.ravel())))

   # Rough estimates of background level and rms:
   sciback = np.sqrt(np.median(var[good]))
   scirms = np.sqrt(np.median(var[good]))
   
   # Report on progress so far:
   if opt.verbose:
      print 'Sci header:', hdr
      print 'Read in sci image:', sci.shape, sci
      print 'Read in var image:', var.shape, var
      print 'Made mask image:', mask.shape, mask
      print 'Useful variance range:', var[good].min(), var[good].max()
      print 'Useful image median level:', sciback
      print 'Useful image median pixel uncertainty:', scirms

   # -------------------------------------------------------------------------
   # Make a first guess at a PSF - a single circularly symmettric Gaussian 
   # defined on same grid as sci image:

   w = np.array([1.0,1.0])           # amplitude at peak
   mu = np.array([[0.0,0.0],[0.0,0.0]])      # centroid position in pixels 
   cov = np.array([[[1.0,0.0],[0.0,1.0]],[[9.0,0.0],[0.0,9.0]]])             # pixels^2, variance matrices
   psf = tractor.GaussianMixturePSF(w,mu,cov)
      
   # -------------------------------------------------------------------------

   # Photometric calibration from PS1 image - Null, because images are 
   # calibrated to return consistent brightness values.
   photocal = tractor.NullPhotoCal()

   # Set up sky to be varied:
   sky = tractor.ConstantSky(0.0)

   # -------------------------------------------------------------------------
   # Make a first guess at a catalog - 4 point sources. Find centre of 
   # field using fitsWCS:

   wcs = lensfinder.PS1WCS(hdr)
   
   x,y,f = NX/2,NY/2, 100*scirms
   e = 5 # pixels
   srcs = [tractor.PointSource(wcs.pixelToPosition(x+e,y),tractor.Flux(f)),
           tractor.PointSource(wcs.pixelToPosition(x-e,y),tractor.Flux(f)),
           tractor.PointSource(wcs.pixelToPosition(x,y+e),tractor.Flux(f)),
           tractor.PointSource(wcs.pixelToPosition(x,y-e),tractor.Flux(f))]

   # -------------------------------------------------------------------------
   
   # Make a tractor Image object out of all this stuff:
   
   image = tractor.Image(data=sci, invvar=invvar,
				 psf=psf, wcs=wcs, sky=sky, photocal=photocal)

   # Start a tractor, and feed it the catalog one src at a time:

   chug = tractor.Tractor([image])
   for src in srcs:
      chug.addSource(src)
   print 'Obtained a total of', len(chug.catalog), 'sources'

   # Freeze all but the PSF, sky and sources:
   for image in chug.images:
      image.freezeParams('photocal', 'wcs', 'psf')

   # Plot:
   plot_state(chug,'initial')
   print chug.getParamNames()

   # Optimize sources with small initial PSF:
   for i in range(5):
      # dlnp2,X,a = chug.optimizeCatalogAtFixedComplexityStep()
      dlnp2,X,a = chug.opt2()
      plot_state(chug,'step-%02d'%i)
      
   # Freeze the sources and thaw the psfs:
   chug.freezeParam('catalog')
   for image in chug.images:
      image.thawParams('psf')
      
   # Optimize everything that is not frozen:
   for i in range(5,10):
      dlnp2,X,a = chug.opt2()
      plot_state(chug,'step-%02d'%i)
#       # Print PSF parameters:
#       print chug.images[0].psf
      # PSF not being optimized correctly
      
   # -------------------------------------------------------------------------
   
   print "Tractor stopping."
      
   return  

# ============================================================================

def plot_state(t,prefix):
   '''
   Make all the plots we need to assess the state of the tractor. Mainly, 
   a multi-panel figure of image, synthetic image and chi, for each image being 
   modelled.
   
   t is a tractor object, containing a list of images.
   '''
   
   px,py = 2,2
   
   for i,image in enumerate(t.images):
      if image.name is None:
         imname = prefix+str(i)
      else:
         imname = image.name
      scale = np.sqrt(np.median(1.0/image.invvar[image.invvar > 0.0]))
 
      ima = dict(interpolation='nearest', origin='lower',
                     vmin=-3.*scale, vmax=20.*scale)
      chia = dict(interpolation='nearest', origin='lower',
                        vmin=-5., vmax=5.)
      psfa = dict(interpolation='nearest', origin='lower')

      plt.figure(figsize=(5*px+1,5*py))
      plt.clf()
      plt.gray

      plt.subplot(py,px,1)
      plt.imshow(image.data, **ima)
      plt.colorbar()

      model = t.getModelImages()[i]
      plt.subplot(py,px,2)
      plt.imshow(model, **ima)
      plt.colorbar()

      chi = t.getChiImage(i)
      plt.subplot(py,px,3)
      plt.imshow(chi, **chia)
      plt.colorbar()

      psfimage = image.psf.getPointSourcePatch(*model.shape).patch
      plt.subplot(py,px,4)
      plt.imshow(psfimage, **psfa)
      plt.colorbar()

      plt.savefig(imname+'.png')   
   
   return

# ============================================================================

if __name__ == '__main__':
   
   ps1tractor()

