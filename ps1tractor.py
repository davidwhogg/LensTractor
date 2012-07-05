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

Bugs
----

 BUG: PSF not being optimized correctly - missing derivatives?

'''

if __name__ == '__main__':
   import matplotlib
   matplotlib.use('Agg')
   # Fonts, latex:
   matplotlib.rc('font',**{'family':'serif', 'serif':['TimesNewRoman'], 'size':18.0})
   matplotlib.rc('text', usetex=True)

import os
import logging
import numpy as np
import pylab as plt
import pyfits

from astrometry.util.file import *
from astrometry.util.util import Tan
from astrometry.util.pyfits_utils import *

import tractor
from tractor import sdss_galaxy as galaxy
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
   
   name = scifile.replace('_sci.fits','')
   
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
      print 'Image name:', name
      print 'Sci header:', hdr
      print 'Read in sci image:', sci.shape, sci
      print 'Read in var image:', var.shape, var
      print 'Made mask image:', mask.shape, mask
      print 'Useful variance range:', var[good].min(), var[good].max()
      print 'Useful image median level:', sciback
      print 'Useful image median pixel uncertainty:', scirms

   # -------------------------------------------------------------------------
   # Make a first guess at a PSF - a single circularly symmetric Gaussian 
   # defined on same grid as sci image:

   w = np.array([1.0])             # amplitude at peak
   mu = np.array([[0.0,0.0]])      # centroid position in pixels 
   FWHM = 4 # pixels
   var = (FWHM/2.35)**2.0
   cov = np.array([[[var,0.0],[0.0,var]]])  # pixels^2, variance matrix
   psf = tractor.GaussianMixturePSF(w,mu,cov)
      
   # # Double Gaussian alternative:
   # w = np.array([1.0,1.0])           # amplitude at peak
   # mu = np.array([[0.0,0.0],[0.0,0.0]])      # centroid position in pixels 
   # FWHM = 7 # pixels
   # var = (FWHM/2.35)**2.0
   # cov = np.array([[[1.0,0.0],[0.0,1.0]],[[var,0.0],[0.0,var]]])  # pixels^2, variance matrices
   # psf = tractor.GaussianMixturePSF(w,mu,cov)
      
   # -------------------------------------------------------------------------

   # Photometric calibration from PS1 image - Null, because we're working in flux
   # units, not calibrated mags; this corresponds to using tractor.Flux() below
   # when creating the Source objects.
   photocal = tractor.NullPhotoCal()

   # Set up sky to be varied:
   sky = tractor.ConstantSky(0.0)

   # -------------------------------------------------------------------------
   # Make a first guess at a catalog - 4 point sources. Find centre of 
   # field using fitsWCS:

   wcs = lensfinder.PS1WCS(hdr)
   
   # Test: 4 point sources:
   #  x,y,f = NX/2,NY/2, 100*scirms
   #  e = 1 # pixels
   #  srcs = [tractor.PointSource(wcs.pixelToPosition(x+e,y),tractor.Flux(f)),
   #          tractor.PointSource(wcs.pixelToPosition(x-e,y),tractor.Flux(f)),
   #          tractor.PointSource(wcs.pixelToPosition(x,y+e),tractor.Flux(f)),
   #          tractor.PointSource(wcs.pixelToPosition(x,y-e),tractor.Flux(f))]
   
   # Source:
   xs,ys, ms = 0.5*NX, 0.5*NY, tractor.Mag(21.0)
   print ms
   sourcepos = wcs.pixelToPosition(xs,ys)
   print sourcepos
   
   pointsource = tractor.PointSource(sourcepos,ms)
   print pointsource
   
   # Lens mass:
   thetaE = lensfinder.EinsteinRadius(1.5) # arcsec
   print thetaE
   gamma = 0.2 # to make quad
   phi   = 0.0 # deg
   xshear = lensfinder.ExternalShear(gamma,phi)
   print xshear
   
   # Lens light:
   x,y = 0.5*NX,0.5*NY
   lenspos = wcs.pixelToPosition(x,y)
   md = tractor.Mag(20.0)
   print md
   re = 1.0  # arcsec
   q = 1.0   # axis ratio
   theta = 0.0 # degrees
   galshape = galaxy.GalaxyShape(re,q,theta)
   print galshape
      
   lensgalaxy = lensfinder.LensGalaxy(lenspos,md,galshape,thetaE,xshear)
   print lensgalaxy
   
   psl = lensfinder.PointSourceLens(lensgalaxy, pointsource)
   print psl
   
   assert False
   
   # srcs = [psl]


   # -------------------------------------------------------------------------
   
   # Make a tractor Image object out of all this stuff:
   
   initimage = tractor.Image(data=sci, invvar=invvar, name=name,
				 psf=psf, wcs=wcs, sky=sky, photocal=photocal)

   # Optimization plan:
   
   Nsteps_optimizing_catalog = 20
   Nsteps_optimizing_PSFs = 0 # To save time.

   # Start a tractor, and feed it the catalog one src at a time:

   chug = tractor.Tractor([initimage])
   for src in srcs:
      chug.addSource(src)
   print 'Obtained a total of', len(chug.catalog), 'sources'

   # Plot:
   plot_state(chug,'progress_initial')

   # Freeze the PSF, sky and photocal, leaving the sources:
   print "DEBUGGING: Before freezing, PSF = ",chug.getImage(0).psf
   for image in chug.getImages():
      image.freezeParams('photocal', 'wcs', 'psf')
   print "DEBUGGING: After freezing, PSF = ",chug.getImage(0).psf
   print "DEBUGGING: pars to be optimized are:",chug.getParamNames()

   # Optimize sources with initial PSF:
   for i in range(Nsteps_optimizing_catalog):
      # dlnp2,X,a = chug.optimizeCatalogAtFixedComplexityStep()
      dlnp2,X,a = chug.optimize()
      plot_state(chug,'progress_optimizing-catalog_step-%02d'%i)
      print "DEBUGGING: pars being optimized are:",chug.getParamNames()
                  
   # Freeze the sources and thaw the psfs:
   print "DEBUGGING: Before thawing, PSF = ",chug.getImage(0).psf
   chug.freezeParam('catalog')
   for image in chug.getImages():
      image.thawParams('psf')
   print "DEBUGGING: After thawing, PSF = ",chug.getImage(0).psf
   print "DEBUGGING: pars to be optimized are:",chug.getParamNames()

   # Optimize everything that is not frozen:
   for i in range(Nsteps_optimizing_catalog,Nsteps_optimizing_catalog+Nsteps_optimizing_PSFs):
      dlnp2,X,a = chug.optimize()
      plot_state(chug,'progress_optimizing-psf_step-%02d'%i)
      print "DEBUGGING: pars being optimized are:",chug.getParamNames()
   
   print "DEBUGGING: After optimizing, PSF = ",chug.getImage(0).psf
   # BUG: PSF not being optimized correctly - missing derivatives?
      
   # -------------------------------------------------------------------------
   
   print "Tractor stopping."
      
   return  

# ============================================================================

def plot_state(t,suffix):
   '''
   Make all the plots we need to assess the state of the tractor. Mainly, 
   a multi-panel figure of image, synthetic image and chi, for each image being 
   modelled.
   
   t is a tractor object, containing a list of images.
   '''
   
   px,py = 2,2
   figprops = dict(figsize=(5*px,5*py), dpi=128)                                          # Figure properties
   adjustprops = dict(\
      left=0.05,\
      bottom=0.05,\
      right=0.95,\
      top=0.95,\
      wspace=0.1,\
      hspace=0.1)
        
   
   for i,image in enumerate(t.images):
      if image.name is None:
         imname = suffix+str(i)
      else:
         imname = image.name
      scale = np.sqrt(np.median(1.0/image.invvar[image.invvar > 0.0]))
 
      ima = dict(interpolation='nearest', origin='lower',
                     vmin=-30.*scale, vmax=3.*scale)
      chia = dict(interpolation='nearest', origin='lower',
                        vmin=-5., vmax=5.)
      psfa = dict(interpolation='nearest', origin='lower')

      fig = plt.figure(**figprops)
      fig.subplots_adjust(**adjustprops)
      plt.clf()
      plt.gray()

      plt.subplot(py,px,1)
      plt.imshow(-image.data, **ima)
      # plt.colorbar()
      tidyup_plot()
      plt.title('Observed image')

      model = t.getModelImages()[i]
      plt.subplot(py,px,2)
      plt.imshow(-model, **ima)
      # plt.colorbar()
      tidyup_plot()
      plt.title('Predicted image')

      chi = t.getChiImage(i)
      plt.subplot(py,px,3)
      plt.imshow(-chi, **chia)
      # plt.colorbar()
      tidyup_plot()
      plt.title('Residuals ($\pm 5\sigma$)')

      psfimage = image.psf.getPointSourcePatch(*model.shape).patch
      plt.subplot(py,px,4)
      plt.imshow(-psfimage, **psfa)
      # plt.colorbar()
      tidyup_plot()
      plt.title('PSF')

      plt.savefig(imname+'_'+suffix+'.png')   
   
   return

def tidyup_plot():
   # Turn off the axis ticks and labels:
   ax = plt.gca()
   ax.xaxis.set_ticks([])
   ax.yaxis.set_ticks([])
   return

# ============================================================================

if __name__ == '__main__':
   
   ps1tractor()

