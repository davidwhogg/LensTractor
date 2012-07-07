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
 - Lens model not being optimized, step sizes and derivatives wrong
 - StepSizes need optimizing for rapid start-up
 - Magic numbers for initial mags, needs automating
 - Memory leak: restrict no. of sampling iterations :-(
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
import time

# from astrometry.util.file import *
from astrometry.util import util
# from astrometry.util.pyfits_utils import *

import tractor
import lenstractor
import emcee

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
     -h --help        Print this message
     -v --verbose     Verbose operation
     -s --sample      Sample the posterior PDF instead of optimizing
     -n --no-plots    Do not plot progress

   INPUTS
     *.fits           Deck of postcard images

   OPTIONAL INPUTS

   OUTPUTS
     stdout                       Useful information
     *.png                        Plots in png format
     
     To be implemented:
       lenstractor_progress.log     Logged output
       lenstractor_results.txt      Model comparison results
       lenstractor_lens.cat         Lens model parameters, including lightcurves
       lenstractor_nebula.cat       Nebula model parameters, including lightcurves

   EXAMPLES

     LensTractor.py -v examples/*.fits

   DEPENDENCIES
     * The Tractor     astrometry.net/svn/trunk/projects/tractor
     * emcee           github.com/danfm/emcee
     * astrometry.net  astrometry.net/svn/trunk/util

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
   # Sampling:
   parser.add_option('-s', '--sample', dest='MCMC', action='count', default=False, help='Sample posterior PDF')
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
   images,total_mags,bands = lenstractor.Deal(scifiles,varfiles,SURVEY='PS1',vb=vb)
   
   # -------------------------------------------------------------------------
   # Initialize the catalog.
   
   # Get rough idea of object position from wcs of first image- works
   # well if all images are the same size and well registered!
   wcs = images[0].wcs
   NX,NY = np.shape(images[0].data)
   
   # m0 = 15.0
   # magnitudes = m0*np.ones(len(bandnames))
   # MAGIC initial magnitude. This should be estimated from
   # the total flux in each (background-subtracted) image...
   
   bandnames = np.unique(bands)
   print bandnames
   magnitudes = np.zeros(len(bandnames))
   for i,bandname in enumerate(bandnames):
       index = np.where(bands == bandname)
       magnitudes[i] = np.median(total_mags[index])
   
   if vb: print "Generic initial SED ",dict(zip(bandnames,magnitudes))
   
   
   # This will be in a 2-model loop eventually:
   
   model = 'nebula'
   
   # model = 'lens'
   
   
   if vb: print "Initializing model: "+model
       
   if model == 'nebula':
   
       # srcs = [InitializeNebula(wcs,bandnames)]

       # Nebula - a flexible galaxy plus four point sources, 
       #   5 sources in total. Start them off with equal magnitudes:
       
       x,y = 0.5*NX,0.5*NY
       fudge = 2
       equalmagnitudes = magnitudes + 2.5*np.log10(5*fudge)
     
       # Composite galaxy:
       galpos = wcs.pixelToPosition(x,y)
       mg = tractor.Mags(order=bandnames, **dict(zip(bandnames,equalmagnitudes)))
       re = 1.0    # arcsec
       q = 1.0     # axis ratio
       theta = 0.0 # degrees
       galshape = tractor.sdss_galaxy.GalaxyShape(re,q,theta)       
       nebulousgalaxy = tractor.sdss_galaxy.CompositeGalaxy(galpos,mg,galshape,mg,galshape)
       if vb: print nebulousgalaxy
       
       # plus 4 point sources, arranged in a small cross:
       e = 3.0 # pixels
       mags = tractor.Mags(order=bandnames, **dict(zip(bandnames,equalmagnitudes)))
       srcs = [nebulousgalaxy,
               tractor.PointSource(wcs.pixelToPosition(x+e,y),mags),
               tractor.PointSource(wcs.pixelToPosition(x-e,y),mags),
               tractor.PointSource(wcs.pixelToPosition(x,y+e),mags),
               tractor.PointSource(wcs.pixelToPosition(x,y-e),mags)]
       
       # Old setup - just 4 point sources:
       # srcs = [tractor.PointSource(wcs.pixelToPosition(x+e,y),mags),
       #         tractor.PointSource(wcs.pixelToPosition(x-e,y),mags),
       #         tractor.PointSource(wcs.pixelToPosition(x,y+e),mags),
       #         tractor.PointSource(wcs.pixelToPosition(x,y-e),mags)]
       

   elif model == 'lens':
          
       # srcs = [InitializeLens(wcs,bandnames)]
       
       # Source to be lensed:
       xs,ys = 0.5*NX, 0.5*NY
       unlensedmagnitudes = magnitudes + 2.5*np.log10(10.0)
       ms = tractor.Mags(order=bandnames, **dict(zip(bandnames,unlensedmagnitudes)))
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
       if vb: print lenspos
       halfmagnitudes = magnitudes + 2.5*np.log10(2.0)
       md = tractor.Mags(order=bandnames, **dict(zip(bandnames,halfmagnitudes)))
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
      print "Model =",srcs
      print " "

   # -------------------------------------------------------------------------

   # Set up logging to theterminal by The Tractor:   
   if opt.verbose:
      lvl = logging.DEBUG
   else:
      lvl = logging.INFO
   logging.basicConfig(level=lvl, format='%(message)s', stream=sys.stdout)

   # -------------------------------------------------------------------------

   # Start a tractor, and let it make a catalog one src at a time:
   chug = tractor.Tractor(images)
   for src in srcs:
      chug.addSource(src)

   # Plot initial state:
   lenstractor.Plot_state(chug,model+'_progress_initial')

   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if not opt.MCMC:
   # Optimize the model parameters:
   
         Nsteps_optimizing_catalog = 20
         Nsteps_optimizing_PSFs = 1    # To save time while developing

         if vb: 
            print "Optimizing model:"
            print "   - no. of iterations to be spent on catalog: ",Nsteps_optimizing_catalog
            print "   - no. of iterations to be spent on PSFs: ",Nsteps_optimizing_PSFs

         # Freeze the PSF, sky and photocal, leaving the sources:
         print "DEBUGGING: Before freezing, PSF = ",chug.getImage(0).psf
         for image in chug.getImages():
            image.freezeParams('photocal', 'wcs', 'psf')
         print "DEBUGGING: After freezing, PSF = ",chug.getImage(0).psf
         print "DEBUGGING: Catalog parameters to be optimized are:",chug.getParamNames()
         
         # THIS LINE FAILS FOR LENS MODEL:
         print "DEBUGGING: Step sizes:",chug.getStepSizes()

         # Optimize sources with initial PSF:
         for i in range(Nsteps_optimizing_catalog):
            dlnp,X,a = chug.optimize()
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
            dlnp,X,a = chug.optimize()
            if not opt.noplots: lenstractor.Plot_state(chug,model+'_progress_optimizing-psf_step-%02d'%i)
            print "DEBUGGING: pars being optimized are:",chug.getParamNames()

         print "DEBUGGING: After optimizing, PSF = ",chug.getImage(0).psf

         # BUG: PSF not being optimized correctly - missing derivatives?

   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   elif opt.MCMC:
   # MCMC sample the model parameters.

         if vb: 
            print "Sampling model parameters with emcee:"
 
         # Freeze the sky and photocal, leaving the PSFs and sources:
         print "DEBUGGING: Before freezing, PSF = ",chug.getImage(0).psf
         for image in chug.getImages():
            image.freezeParams('photocal', 'wcs')
            
         # Get the thawed parameters:
         p0 = np.array(chug.getParams())
         print 'Tractor parameters:'
         for i,parname in enumerate(chug.getParamNames()):
               print '  ', parname, '=', p0[i]
         ndim = len(p0)
         print 'Number of parameter space dimensions: ',ndim

         # Make an emcee sampler that uses our tractor to compute its logprob:
         nw = 8*ndim
         sampler = emcee.EnsembleSampler(nw, ndim, chug, live_dangerously=True)

         # Start the walkers off near the initialisation point - 
         # We need it to be ~1 pixel in position, and not too much
         # flux restrction... 
         
         if model=='lens':
            # The following gets us 0.2" in dec:
            psteps = np.zeros_like(p0) + 0.00004
            # This could be optimized, to allow more initial freedom in eg flux.
         
         elif model=='nebula':
            # Good first guess should be some fraction of the optimization step sizes:
            psteps = 0.2*np.array(chug.getStepSizes())
         
         print "Initial size (in each dimension) of sample ball = ",psteps
         
         pp = emcee.EnsembleSampler.sampleBall(p0, psteps, nw)
         rstate = None
         lnp = None
         
         # Take a few steps - memory leaks fast! (~10Mb per sec)
         for step in range(1, 4):
               
               print 'Run MCMC step set:', step
               t0 = tractor.Time()
               pp,lnp,rstate = sampler.run_mcmc(pp, 5, lnprob0=lnp, rstate0=rstate)
               
               print 'Mean acceptance fraction after', sampler.iterations, 'iterations =',np.mean(sampler.acceptance_fraction)
               t_mcmc = (tractor.Time() - t0)
               print 'Runtime:', t_mcmc
               
               # Find the current posterior means:
               pbar = np.mean(pp,axis=0)
               print "Mean parameters: ",pbar,np.mean(lnp)
               
               # Find the current best sample:
               maxlnp = np.max(lnp)
               best = np.where(lnp == maxlnp)
               pbest = np.ravel(pp[best,:])
               print "Best parameters: ",pbest,maxlnp
                             
               if not opt.noplots: 
                  chug.setParams(pbest)
                  lenstractor.Plot_state(chug,model+'_progress_sampling_step-%02d'%step)

         print 'Best lnprob:', np.max(lnp)
         # print 'dlnprobs:', ', '.join(['%.1f' % d for d in lnp - np.max(lnp)])
         # print 'MCMC took', t_mcmc, 'sec'

         # Take the last best sample and call it a result:
         chug.setParams(pbest)

   # -------------------------------------------------------------------------
   
   lenstractor.Plot_state(chug,model+'_progress_zcomplete')
   
   print "Tractor stopping."
      
   return  

# ============================================================================

if __name__ == '__main__':
   
   main()

