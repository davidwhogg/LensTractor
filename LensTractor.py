'''
This file is part of the lenstractor project.
Copyright 2012 David W. Hogg (NYU) and Phil Marshall (Oxford).

Description
-----------

Run the Tractor on a deck of single object cutout images.
Read in an image and its weight map, guess the PSF, put an object at the
centre image of the field and then optimize the catalog and PSF.

Example use
-----------

 python LensTractor.py -v examples/H1413+117*.fits

Bugs
----

 - PSF not being optimized correctly - missing derivatives?
 - Header PSF FWHM sometimes NaN, hard to recover from
 - lens model not being optimized, step sizes and derivatives wrong

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
import pyfits

# from astrometry.util.file import *
from astrometry.util import util
# from astrometry.util.pyfits_utils import *

import tractor
import lenstractor

# ============================================================================

def main():

   """
   NAME
     LensTractor.py

   PURPOSE
     Run the Tractor on a deck of single object cutout images.
     Read in an image and its weight map, guess the PSF, put an object at the
     centre image of the field and then optimize the catalog and PSF.

   COMMENTS

   FLAGS
     -h            Print this message [0]
     -v            Verbose operation [0]

   INPUTS
     *.fits        Deck of postcard images

   OPTIONAL INPUTS
     --no-plots    Do not plot progress [0]

   OUTPUTS
     stdout                       Useful information
     *.png                        Plots in png format
     lenstractor_progress.log     Logged output
     lenstractor_results.txt      Model comparison results
     lenstractor_lens.cat         Lens model parameters, including lightcurves
     lenstractor_nebula.cat         Nebula model parameters, including lightcurves

   EXAMPLES

     LensTractor.py -v examples/*.fits

   BUGS
     - PSFs not being optimized correctly - missing derivatives?

   HISTORY
     2012-07-06       First predicted Lens images Marshall/Hogg (Oxford/NYU)
   """

   # --------------------------------------------------------------------

   from optparse import OptionParser
   import sys

   # Set available options:
   parser = OptionParser(usage=('%prog *.fits'))
   # Verbosity:
   parser.add_option('-v', '--verbose', dest='verbose', action='count', default=False, help='Make more verbose')
   vb = True  # for usual outputs.
   # Plotting:
   parser.add_option('-n', '--no-plots', dest='noplots', action='count', default=False, help='Skip plotting')
   
   # Read in options and arguments - note only sci and wht images are supplied:
   opt,args = parser.parse_args()
   
   if len(args) < 2:
      parser.print_help()
      sys.exit(-1)
    
   # The rest of the command line is assumed to be a list of files:
   inputfiles = args
 
   # -------------------------------------------------------------------------
   
   # Organise the deck of inputfiles into scifiles and varfiles:
   scifiles,varfiles = lenstractor.Riffle(inputfiles,vb=vb)
   
   # Read into Tractor Image objects, and see what filters we have:   
   images,bandnames,epochs = lenstractor.Deal(scifiles,varfiles,SURVEY='PS1',vb=vb)
   
   # -------------------------------------------------------------------------
   # Initialize the catalog.
   
   # Get rough idea of object position from wcs first image.
   wcs = images[0].wcs
   NX,NY = np.shape(images[0].data)
   
   # This will be in a 2-model loop eventually:
   # model = 'nebula'
   model = 'lens'
   
   if vb: print "Initializing model: "+model
       
   if model == 'nebula':
   
       # srcs = [InitializeNebula(wcs,bandnames)]

       # Placeholder: 4 point sources:
       x,y = NX/2.0,NY/2.0
       e = 3.0 # pixels
       magnitudes = 16.0*np.ones(len(bandnames))
       if vb: print "Initial SED ",dict(zip(bandnames,magnitudes))
       
       mags = tractor.Mags(order=bandnames, **dict(zip(bandnames,magnitudes)))
       srcs = [tractor.PointSource(wcs.pixelToPosition(x+e,y),mags),
               tractor.PointSource(wcs.pixelToPosition(x-e,y),mags),
               tractor.PointSource(wcs.pixelToPosition(x,y+e),mags),
               tractor.PointSource(wcs.pixelToPosition(x,y-e),mags)]
       
 
   elif model == 'lens':
          
       # srcs = [InitializeLens(wcs,bandnames)]
       
       # Source to be lensed:
       xs,ys = 0.5*NX, 0.5*NY
       magnitudes = 17.0*np.ones(len(bandnames))
       if vb: print "Initial source SED ",dict(zip(bandnames,magnitudes))
       ms = tractor.Mags(order=bandnames, **dict(zip(bandnames,magnitudes)))
       if vb: print ms
       sourcepos = wcs.pixelToPosition(xs,ys)
       if vb: print sourcepos

       pointsource = tractor.PointSource(sourcepos,ms)
       if vb: print pointsource

       # Lens mass:
       thetaE = lenstractor.EinsteinRadius(0.75) # arcsec
       if vb: print thetaE
       gamma = 0.2 # to make quad
       phi   = 0.0 # deg
       xshear = lenstractor.ExternalShear(gamma,phi)
       if vb: print xshear

       # Lens light:
       x,y = 0.5*NX,0.5*NY
       lenspos = wcs.pixelToPosition(x,y)
       magnitudes = 17.0*np.ones(len(bandnames))
       if vb: print "Initial lens SED ",dict(zip(bandnames,magnitudes))
       md = tractor.Mags(order=bandnames, **dict(zip(bandnames,magnitudes)))
       if vb: print md
       re = 1.0  # arcsec
       q = 1.0   # axis ratio
       theta = 0.0 # degrees
       galshape = tractor.sdss_galaxy.GalaxyShape(re,q,theta)
       if vb: print galshape

       lensgalaxy = lenstractor.LensGalaxy(lenspos,md,galshape,thetaE,xshear)
       if vb: print lensgalaxy

       srcs = [lenstractor.PointSourceLens(lensgalaxy, pointsource)]

   if vb: 
      print srcs
      print " "

   # -------------------------------------------------------------------------

   # Set up logging to theterminal by The Tractor:   
   if opt.verbose:
      lvl = logging.DEBUG
   else:
      lvl = logging.INFO
   logging.basicConfig(level=lvl, format='%(message)s', stream=sys.stdout)

   # -------------------------------------------------------------------------
   # Optimize the model parameters.
   
   Nsteps_optimizing_catalog = 5
   Nsteps_optimizing_PSFs = 2    # To save time while developing
   
   if vb: 
      print "Optimizing model:"
      print "   - no. of iterations to be spent on catalog: ",Nsteps_optimizing_catalog
      print "   - no. of iterations to be spent on PSFs: ",Nsteps_optimizing_PSFs

   # Start a tractor, and let it make a catalog one src at a time:
   chug = tractor.Tractor(images)
   for src in srcs:
      chug.addSource(src)

   # Plot initial state:
   if not opt.noplots: lenstractor.Plot_state(chug,model+'_progress_initial')

   # Freeze the PSF, sky and photocal, leaving the sources:
   print "DEBUGGING: Before freezing, PSF = ",chug.getImage(0).psf
   for image in chug.getImages():
      image.freezeParams('photocal', 'wcs', 'psf')
   print "DEBUGGING: After freezing, PSF = ",chug.getImage(0).psf
   print "DEBUGGING: Catalog parameters to be optimized are:",chug.getParamNames()
   print "DEBUGGING: Step sizes:",chug.getStepSizes()

   # Optimize sources with initial PSF:
   for i in range(Nsteps_optimizing_catalog):
      # dlnp2,X,a = chug.optimizeCatalogAtFixedComplexityStep()
      dlnp2,X,a = chug.optimize()
      if not opt.noplots: lenstractor.Plot_state(chug,model+'_progress_optimizing-catalog_step-%02d'%i)
      print "DEBUGGING: pars being optimized are:",chug.getParamNames()
                  
   # BUG: lens not being optimized correctly - missing derivatives?
      
   # Freeze the sources and thaw the psfs:
   print "DEBUGGING: Before thawing, PSF = ",chug.getImage(0).psf
   chug.freezeParam('catalog')
   for image in chug.getImages():
      image.thawParams('psf')
   print "DEBUGGING: After thawing, PSF = ",chug.getImage(0).psf
   print "DEBUGGING: PSF parameters to be optimized are:",chug.getParamNames()
   print "DEBUGGING: Step sizes:",chug.getStepSizes()

   # Optimize everything that is not frozen:
   for i in range(Nsteps_optimizing_catalog,Nsteps_optimizing_catalog+Nsteps_optimizing_PSFs):
      dlnp2,X,a = chug.optimize()
      if not opt.noplots: lenstractor.Plot_state(chug,model+'_progress_optimizing-psf_step-%02d'%i)
      print "DEBUGGING: pars being optimized are:",chug.getParamNames()
   
   print "DEBUGGING: After optimizing, PSF = ",chug.getImage(0).psf
   
   # BUG: PSF not being optimized correctly - missing derivatives?
      
   # -------------------------------------------------------------------------
   
   lenstractor.Plot_state(chug,model+'_progress_complete')
   
   print "Tractor stopping."
      
   return  

# ============================================================================

if __name__ == '__main__':
   
   main()

