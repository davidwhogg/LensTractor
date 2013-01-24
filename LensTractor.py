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

 python LensTractor.py -n examples/H1413+117*.fits

Bugs
----
 - Lens initialisation, esp source positions, needs careful attention
 - StepSizes need optimizing for lens model, esp source position

 - Point source mags are not variable
 
 - PSF not being optimized correctly - missing derivatives?
 - Header PSF FWHM sometimes NaN, no recovery from this yet
 - Memory leak: restrict no. of sampling iterations :-(
 - PhotoCal may need optimizing if zpts are untrustworthy!
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
     The idea is to identify good lens candidates by principled model 
     selection: two well-defined models competing against each other, given 
     multi-epoch imaging data. The Nebula model (1 extended source, plus
     N=1,2,3 or 4 point sources, with sub-models denoted by "NebulaN") is very
     flexible, so should be better at fitting the data in general than the
     Lens model (1 extended source, plus 1 background point source). However,
     when the Lens provides a good fit, it does so at lower cost (fewer
     parameters), so should win by Bayesian information criteria (we use BIC
     as a cheap proxy for evidence ratio).
     
     The default workflow is as follows:
     
       Is the target a lens?
       
       * Try Nebula1
       
       * Try Nebula2
       
           if Nebula1 beats Nebula2: 
             Return NO
           else:
             Nebula = Nebula2
       
       * Try Nebula4
       
           if Nebula4 beats Nebula2: 
             Nebula = Nebula4

                      
       * Try Lens (via Nebula)
           if Lens beats Nebula: 
             Return YES
           else:
             Return NO

      Initialisation of Lens via Nebula depends on the point source 
      multiplicity of the final Nebula model: what we do with three point
      image positions will be different from what we do with 4 point image
      positions, particularly with regard to the deflector centroid. 
      
      Open questions:
      
      Does it make sense to dogmatically associate the extended object with
      the deflector?
         YES: detection of a deflector galaxy will be needed for a convincing
         candidate anyway.
         NO: using the extended object to model a high S/N merging image
         system should not be punished
      
      How are we going to interpret the point image positions if we do not
      have an estimated deflector position?
      

   FLAGS
     -h --help        Print this message
     -v --verbose     Verbose operation
     -s --sample      Sample the posterior PDF instead of optimizing
     -x --no-plots    Do not plot progress
     -l --lens        Only fit lens model, initialized from scratch
     -n --nebula      Only fit nebula model, initialized from scratch

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

     python LensTractor.py -x examples/H1413+117_10x10arcsec_55*fits > examples/H1413+117_10x10arcsec_lenstractor.log
   
   DEPENDENCIES
     * The Tractor     astrometry.net/svn/trunk/projects/tractor
     * emcee           github.com/danfm/emcee
     * astrometry.net  astrometry.net/svn/trunk/util

   BUGS
     - PSFs not being optimized correctly - permafrost?

   HISTORY
     2012-07-06       First predicted Lens images Marshall/Hogg (Oxford/NYU)
     2013-01-23       Sequential Nebula 1,2,4 scheme implemented
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
   parser.add_option('-x', '--no-plots', dest='noplots', action='count', default=False, help='Skip plotting')
   # Lens model only:
   parser.add_option('-l', '--lens', dest='lens', action='count', default=False, help='Fit lens model')
   # Nebula model only:
   parser.add_option('-n', '--nebula', dest='nebula', action='count', default=False, help='Fit nebula model')
   
   # Read in options and arguments - note only sci and wht images are supplied:
   opt,args = parser.parse_args()
   
   if len(args) < 2:
      parser.print_help()
      sys.exit(-1)
    
   # The rest of the command line is assumed to be a list of files:
   inputfiles = args
   
   # Workflow:
   if opt.lens:
      models = ['lens']
   elif opt.nebula:
      models = ['nebula1','nebula2','nebula4',]
   else:
      models = ['nebula1','nebula2','nebula4','lens']
   BIC = dict(zip(models,np.zeros(len(models))))
   # NB. default operation is to fit both and compare.
   
   if vb: 
      print "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
      print "                               LensTractor "
      print "    Fitting",models," models to a deck of FITS postcards"
      print "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
 
   # -------------------------------------------------------------------------
   
   # Organise the deck of inputfiles into scifiles and varfiles:
   scifiles,varfiles = lenstractor.Riffle(inputfiles,vb=vb)
   
   # Read into Tractor Image objects, and see what filters we have:   
   images,total_mags,bands = lenstractor.Deal(scifiles,varfiles,SURVEY='PS1',vb=vb)
   
   # -------------------------------------------------------------------------
   # Generic items needed to initialize the Tractor's catalog.
   
   # Get rough idea of object position from wcs of first image- works
   # well if all images are the same size and well registered!
   wcs = images[0].wcs
   NX,NY = np.shape(images[0].data)
   if vb: print "Generic initial position ",NX,NY,"(pixels)"
   
   # m0 = 15.0
   # magnitudes = m0*np.ones(len(bandnames))
   # MAGIC initial magnitude. This should be estimated from
   # the total flux in each (background-subtracted) image...
   
   bandnames = np.unique(bands)
   magnitudes = np.zeros(len(bandnames))
   for i,bandname in enumerate(bandnames):
       index = np.where(bands == bandname)
       magnitudes[i] = np.median(total_mags[index])
   if vb: print "Generic initial SED ",dict(zip(bandnames,magnitudes))
   
   # -------------------------------------------------------------------------
   
   # Loop over models:
      
   for model in models: 
      
       if vb: 
           print "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
           print "Initializing model: "+model

       # Figure out what type of model this is:
       modeltype =  model[0:6]
       
       # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       
       if modeltype == 'nebula':

           # Nebula - a flexible galaxy plus K point sources, 

           # How many point sources?
           K = int(model[6:7])
           
           # The first nebula model we run has 1 point source and one 
           # galaxy, initialised sensibly but randomly. 
           # All subsequent models just have extra random point sources.

           # Start both sources off with equal magnitudes:
           x,y = 0.5*NX,0.5*NY
           fudge = 2
           equalmagnitudes = magnitudes + 2.5*np.log10(5*fudge)
           mags = tractor.Mags(order=bandnames, **dict(zip(bandnames,equalmagnitudes)))

           # Add an exponential galaxy:
           galpos = wcs.pixelToPosition(x,y)
           mg = tractor.Mags(order=bandnames, **dict(zip(bandnames,equalmagnitudes)))
           re = 0.5    # arcsec
           q = 1.0     # axis ratio
           theta = 0.0 # degrees
           galshape = tractor.sdss_galaxy.GalaxyShape(re,q,theta)       
           nebulousgalaxy = tractor.sdss_galaxy.ExpGalaxy(galpos,mags.copy(),galshape)
           if vb: print nebulousgalaxy
           srcs = [nebulousgalaxy]

           for i in range(K):
               # Add a point source with random position near nebula centre:
               e = 2.0 # pixels
               x,y = e*np.random.randn(2)
               star = tractor.PointSource(wcs.pixelToPosition(x+e,y+e),mags.copy())
               if vb: print star
               srcs.append(star)

       # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

       # Original Nebula4 initialisation went as follows:
       #   tractor.PointSource(wcs.pixelToPosition(x+e,y),mags.copy()),
       #   tractor.PointSource(wcs.pixelToPosition(x-e,y),mags.copy()),
       #   tractor.PointSource(wcs.pixelToPosition(x,y+e),mags.copy()),
       #   tractor.PointSource(wcs.pixelToPosition(x,y-e),mags.copy())]

       # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       
       elif modeltype == 'lens':

           # Source to be lensed:
           xs,ys = 0.5*NX, 0.5*NY
           unlensedmagnitudes = magnitudes + 2.5*np.log10(40.0)
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

       # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       
       if vb: 
          print "Model =",srcs
          print " "

       # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

       # Set up logging to the terminal by The Tractor:   
       if opt.verbose:
          lvl = logging.DEBUG
       else:
          lvl = logging.INFO
       logging.basicConfig(level=lvl, format='%(message)s', stream=sys.stdout)

       # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

       # Start a tractor, and let it make a catalog one src at a time.
       # Pass in a copy of the image list, so that PSF etc are 
       # initialised correctly for each model. 
       chug = tractor.Tractor(list(images))
       for src in srcs:
          chug.addSource(src)

       # Plot initial state:
       lenstractor.Plot_state(chug,model+'_progress_initial')

       # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

       if not opt.MCMC:
       # Optimize the model parameters:

             if modeltype=='nebula':
                Nrounds = 1 # 2
                Nsteps_optimizing_catalog = 10 # 20
                Nsteps_optimizing_PSFs = 1 # 10
             elif modeltype=='lens':
                Nrounds = 3
                Nsteps_optimizing_catalog = 7
                Nsteps_optimizing_PSFs = 3

             if vb: 
                print "Optimizing model:"
                print "   - no. of iterations per round to be spent on catalog: ",Nsteps_optimizing_catalog
                print "   - no. of iterations per round to be spent on PSFs: ",Nsteps_optimizing_PSFs
                print "   - no. of rounds: ",Nrounds

             k = 0
             for round in range(Nrounds):
                   
                   print "Fitting "+model+": seconds out, round",round
            
                   # Freeze the PSF, sky and photocal, leaving the sources:
                   chug.thawParam('catalog')
                   for image in chug.getImages():
                     image.thawParams('sky')
                     image.freezeParams('photocal', 'wcs', 'psf')
                   print "Fitting "+model+": Catalog parameters to be optimized are:",chug.getParamNames()
                   print "Fitting "+model+": Initial values are:",chug.getParams()
                   print "Fitting "+model+": Step sizes:",chug.getStepSizes()

                   # Optimize sources with initial PSF:
                   for i in range(Nsteps_optimizing_catalog):
                      dlnp,X,a = chug.optimize(damp=3)
                      if not opt.noplots: lenstractor.Plot_state(chug,model+'_progress_optimizing_step-%02d_catalog'%k)
                      print "Fitting "+model+": at step",k,"parameter values are:",chug.getParams()
                      k += 1

                   # Freeze the sources and thaw the psfs:
                   chug.freezeParam('catalog')
                   for image in chug.getImages():
                      image.thawParams('psf')
                      image.freezeParams('photocal', 'wcs', 'sky')
                   print "Fitting PSF: After thawing, zeroth PSF = ",chug.getImage(0).psf
                   print "Fitting PSF: PSF parameters to be optimized are:",chug.getParamNames()
                   print "Fitting PSF: Initial values are:",chug.getParams()
                   print "Fitting PSF: Step sizes:",chug.getStepSizes()

                   # Optimize everything that is not frozen:
                   for i in range(Nsteps_optimizing_PSFs):
                      dlnp,X,a = chug.optimize()
                      if not opt.noplots: lenstractor.Plot_state(chug,model+'_progress_optimizing_step-%02d_catalog'%k)
                      print "Fitting PSF: at step",k,"parameter values are:",chug.getParams()
                      k += 1
                   print "Fitting PSF: After optimizing, zeroth PSF = ",chug.getImage(0).psf

             # BUG: PSF not being optimized correctly - missing derivatives?

             lenstractor.Plot_state(chug,model+'_progress_optimizing_zcomplete')
             

       # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

       elif opt.MCMC:
       # MCMC sample the model parameters.

             if vb: 
                print "Sampling model parameters with emcee:"

             # Freeze the sky and photocal, leaving the PSFs and sources:
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

             lenstractor.Plot_state(chug,model+'_progress_sampling_zcomplete')

       if vb: 
           print "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
       
       # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       
       # Collect statistics about this model's fit:
       
       chisq = -2.0*chug.getLogLikelihood()
       chug.thawParam('catalog')
       for image in chug.getImages():
         image.thawParams('sky', 'psf')
         image.freezeParams('photocal', 'wcs')
       K = len(chug.getParams())
       N = chug.getNdata()
       BIC[model] = chisq + K*np.log(1.0*N)
       print "Fitting "+model+": chisq, K, N, BIC =",chisq,K,N,BIC[model]
       
   # -------------------------------------------------------------------------
   
   if len(models) == 2:
   # Compare models and report:
       print "BIC = ",BIC
       print "Fitting result: Bayes factor in favour of nebula is exp[",-0.5*(BIC['nebula'] - BIC['lens']),"]"
       print "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"

   # -------------------------------------------------------------------------
   
   print "Tractor stopping."
      
   return  

# ============================================================================

if __name__ == '__main__':
   
   main()

