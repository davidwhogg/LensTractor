'''
This file is part of the lenstractor project.
Copyright 2012 David W. Hogg (NYU) and Phil Marshall (Oxford).
'''
# ============================================================================

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


from astrometry.util import util

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
     The idea is to identify good lens candidates by principled model 
     selection: two well-defined models competing against each other, given 
     multi-epoch imaging data. The Nebula model (1 extended source, plus
     N=2, 3 or 4 point sources, with sub-models denoted by "NebulaN") is very
     flexible, so should be better at fitting the data in general than the
     Lens model (1 extended source, plus 1 background point source). However,
     when the Lens provides a good fit, it does so at lower cost (fewer
     parameters), so should win by Bayesian information criteria (we use BIC
     as a cheap proxy for evidence ratio).
     
     The workflow we probably want to aim for is something like the following:
     
       * Fit PSF images with PSF models; fix PSFs
       
       * Try Nebula2
       
       * Try Nebula4 
           if Nebula2 beats Nebula4: 
             Nebula = Nebula2
           else:
             Nebula = Nebula4
                      
       * Try Lens (initialised with Nebula)
           if Lens beats Nebula: 
             Classification = 'Lens'
             Return YES
           else:
             Classification = 'Nebula'
             Return NO

      Initialisation of Lens via Nebula could be tricky - there is some 
      parsing to be done, and decisions to be made... In practice we may end
      up working just with the Nebula output, which should be at least 
      easier to interpret than a SExtractor catalog, for example.
      
      Open questions:
      
      Does it make sense to dogmatically associate the extended object with
      the deflector?
         YES: detection of a deflector galaxy will be needed for a convincing
         candidate anyway.
         NO: using the extended object to model a high S/N merging image
         system should not be punished
      To be investigated.
      
      How are we going to interpret the point image positions if we do not
      have an estimated deflector position?
      

   OPTIONS
     -h --help        Print this message
     -v --verbose     Verbose operation
     -s --sample      Sample the posterior PDF instead of optimizing
     -x --no-plots    Do not plot progress
     -l --lens        Only fit lens model, initialized from scratch

   INPUTS
     *.fits           Deck of postcard images

   OPTIONAL INPUTS
     -n --nebula                  K    Only fit NebulaK model, initialized from scratch
     --optimization-rounds        Nr   Number of rounds of optimization [2]
     --optimization-steps-catalog Nc   Number of steps per round spent
                                        optimizing source catalog [10]
     --optimization-steps-psf     Np   Number of steps per round spent
                                        optimizing PSF catalog [2]
     -o --output             outfile   Name of output catalog filename       

   OUTPUTS
     stdout                       Useful information
     *.png                        Plots in png format
     
     To be implemented:
       lenstractor_progress.log     Logged output
       lenstractor_results.txt      Model comparison results
       lenstractor_lens.cat         Lens model parameters, including lightcurves
       lenstractor_nebula.cat       Nebula model parameters, including lightcurves

   EXAMPLES

     python LensTractor.py -n 4 \
       -o examples/ps1/H1413+117_10x10arcsec_Nebula4.cat \
          examples/ps1/H1413+117_10x10arcsec_55*fits > \
          examples/ps1/H1413+117_10x10arcsec_Nebula4.log
       
     python LensTractor.py -n 2 \
       -o examples/sdss/0951+2635/0951+2635_Nebula2.cat \
          examples/sdss/0951+2635/*sci.fits > \
          examples/sdss/0951+2635/0951+2635_Nebula2.log

       
   
   DEPENDENCIES
     * The Tractor     astrometry.net/svn/trunk/projects/tractor
     * emcee           github.com/danfm/emcee
     * astrometry.net  astrometry.net/svn/trunk/util

   BUGS
     - SDSS examples show bad WCS treatment...
     - Possible problems with image alignment
     - Memory leak: restrict no. of sampling iterations :-(
     - Header PSF FWHM sometimes NaN, no recovery from this yet
     
   FEATURE REQUESTS  
     - Lens initialisation, esp source positions, needs careful attention
     - StepSizes need optimizing for lens model, esp source position 
     - Point source mags are not variable
     - PSF not being optimized correctly - missing derivatives?
     - PhotoCal may need optimizing if zpts are untrustworthy!

   HISTORY
     2012-07-06       First predicted Lens images Marshall/Hogg (Oxford/NYU)
     2013-08-         Adapted for KIDS Buddelmeier (Kapteyn)
     2014-04-         Refactoring for easier experimentation Marshall/Agnello (KIPAC/UCSB)
   """

   # --------------------------------------------------------------------

   from argparse import ArgumentParser
   import sys

   # Set available options:
   parser = ArgumentParser()
   # List of files:
   parser.add_argument('inputfiles', metavar='N', nargs='+')
   # Verbosity:
   parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', default=False, help='Make more verbose')
   # Sampling:
   parser.add_argument('-s', '--sample', dest='MCMC', action='store_true', default=False, help='Sample posterior PDF')
   # Plotting:
   parser.add_argument('-x', '--no-plots', dest='noplots', action='store_true', default=False, help='Skip plotting')
   # Lens model only:
   parser.add_argument('-l', '--lens', dest='lens', action='store_true', default=False, help='Fit Lens model')
   # Nebula model only:
   parser.add_argument('-n', '--nebula', dest='K', type=int, default=0, help='Fit NebulaK model, provide K')
   # Output filename:
   parser.add_argument('-o', '--output', dest='outfile', type=str, default='lenstractor.cat', help='Output catalog filename')
   # Optimization workflow:
   parser.add_argument('--optimization-rounds', dest='Nr', type=int, default=3, help='No. of optimization rounds')
   parser.add_argument('--optimization-steps-catalog', dest='Nc', type=int, default=5, help='No. of optimization steps spent on source catalog')
   parser.add_argument('--optimization-steps-psf', dest='Np', type=int, default=0, help='No. of optimization steps spent on PSFs')
   parser.add_argument('--survey', dest='survey', type=str, default="PS1", help="Survey (SDSS, PS1 or KIDS)")

   # Read in options and arguments - note only sci and wht images are supplied:
   args = parser.parse_args()
   
      
   if (len(args.inputfiles) < 2):
      # parser.print_help()
      print main.__doc__  # Whoah! What does this do?! Some sort of magic.
      sys.exit(-1)
   
   vb = args.verbose
   
   # Workflow:
   if args.lens:
      models = ['Lens']
   elif args.K > 0:
      models = ['Nebula'+str(args.K)]
      # NB. Global default is Nebula1!
   else:
      models = ['Nebula2','Nebula4','Lens']
         
   BIC = dict(zip(models,np.zeros(len(models))))

   # Package up settings:
   opt_settings = {'Nr':args.Nr, 'Nc':args.Nc, 'Np':args.Np}
   # Magic sampling numbers!
   mcmc_settings = {'nwp':20, 'ns':20, 'nss':100}

   if vb: 
      print "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
      print "                               LensTractor "
      print "    Fitting",models," models to a deck of FITS postcards"
      print "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
 
   # -------------------------------------------------------------------------
   
   # Organise the deck of inputfiles into scifiles and varfiles:
   scifiles,varfiles = lenstractor.Riffle(args.inputfiles,vb=vb)
   
   # Read into Tractor Image objects, and see what filters we have:   
   images,total_mags,bands = lenstractor.Deal(scifiles,varfiles,SURVEY=args.survey,vb=vb)
   
   # Package up:
   dataset = list(images)
   
   # -------------------------------------------------------------------------
   # Generic items needed to initialize the Tractor's catalog.
   
   # Get rough idea of object position from wcs of first image- works
   # well if all images are the same size and well registered!
   wcs = images[0].wcs # HACK! Need to consider different WCS in different images...
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
   
   # Loop over models, initialising and fitting.
      
   for thismodel in models: 
      
       if vb: 
           print "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
           print "Initializing model: "+thismodel       
       
       # Figure out what type of model this is:
       modeltype =  thismodel[0:6]
              
       model = lenstractor.Model(thismodel)

       #   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   - -
              
       if modeltype == 'Nebula':

           # Nebula - a flexible galaxy plus K point sources. 

           # How many point sources?
           K = int(thismodel[6:7])
           
           # The first Nebula model we run has 1 point source and one 
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
           model.srcs.append(nebulousgalaxy)

           for i in range(K):
               # Add a point source with random position near nebula centre:
#               e = 2.0 # pixels
               e = 1.0 # pixels
               dx,dy = e*np.random.randn(2)
               # dx,dy = 2.0*np.random.randn(2)
               # Start in cross formation (for testing):
               if i == 0: dx,dy = -3,0
               if i == 1: dx,dy =  0,3
               if i == 2: dx,dy =  3,0
               if i == 3: dx,dy =  0,-3
               star = tractor.PointSource(wcs.pixelToPosition(x+dx,y+dy),mags.copy())
               if vb: print star
               model.srcs.append(star)
           
       #   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   - -
       
       # Original Nebula4 initialisation went as follows:
       #   tractor.PointSource(wcs.pixelToPosition(x+e,y),mags.copy()),
       #   tractor.PointSource(wcs.pixelToPosition(x-e,y),mags.copy()),
       #   tractor.PointSource(wcs.pixelToPosition(x,y+e),mags.copy()),
       #   tractor.PointSource(wcs.pixelToPosition(x,y-e),mags.copy())]

       #   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   - -
              
       elif modeltype == 'Lens':

           # If nebula has been run, use the best nebula model (by BIC)
           # to initialise the lens model. If it hasn't, do something
           # sensible.
           
           
           # Source to be lensed:
           xs,ys = 0.5*NX, 0.5*NY
           # Tiny random offset (pixels):
           e = 0.0
           dx,dy = e*np.random.randn(2)
           xs,ys = xs+dx,ys+dy
           unlensedmagnitudes = magnitudes + 2.5*np.log10(40.0)
           ms = tractor.Mags(order=bandnames, **dict(zip(bandnames,unlensedmagnitudes)))
           if vb: print ms
           sourcepos = wcs.pixelToPosition(xs,ys)
           if vb: print sourcepos

           pointsource = tractor.PointSource(sourcepos,ms)
           if vb: print pointsource

           # Lens mass:
#           thetaE = lenstractor.EinsteinRadius(0.75) # arcsec
           thetaE = lenstractor.EinsteinRadius(0.2) # arcsec

           if vb: print thetaE
           gamma = 0.2 # to make quad
           phi   = 0.0 # deg
           xshear = lenstractor.ExternalShear(gamma,phi)
           if vb: print xshear

           # Lens light:
           x,y = 0.5*NX,0.5*NY
           lenspos = wcs.pixelToPosition(x,y)
           if vb: print lenspos
           lensmagnitudes = magnitudes + 2.5*np.log10(10.0)
           md = tractor.Mags(order=bandnames, **dict(zip(bandnames,lensmagnitudes)))
           if vb: print md
           re = 0.5  # arcsec
           q = 0.8   # axis ratio
           theta = 90.0 # degrees
           galshape = tractor.sdss_galaxy.GalaxyShape(re,q,theta)
           if vb: print galshape

           lensgalaxy = lenstractor.LensGalaxy(lenspos,md,galshape,thetaE,xshear)
           if vb: print lensgalaxy

           model.srcs.append(lenstractor.PointSourceLens(lensgalaxy, pointsource))


       #   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   - -
       
       if vb: 
           print "Initialization complete."
           print "Model =",model
           print " "

       # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

       # Set up logging to the terminal by The Tractor:   
       if vb:
          # lvl = logging.DEBUG
          lvl = logging.INFO
       else:
          # lvl = logging.INFO
          lvl = logging.ERROR
       logging.basicConfig(level=lvl, format='%(message)s', stream=sys.stdout)

       # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

       # Start a lenstractor, which will make a catalog one src at a time.
       # Pass in a copy of the image list, so that the PSF etc are 
       # initialised correctly for each model. 
       
       LT = lenstractor.LensTractor(dataset,model,args.survey,vb=vb,noplots=args.noplots)

       # Plot initial state:
       lenstractor.Plot_state(
           LT.chug,
           LT.model.name+'_progress_initial',
           SURVEY=args.survey)

       # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

       if not args.MCMC:

           LT.drive(by='optimizing',using=opt_settings)

       else:

           LT.drive(by='sampling',using=mcmc_settings)

       # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       
       if vb:
           print "Fit complete."
           print "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
       
       # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              
       # Compute BIC for this fit:
       BIC[thismodel] = LT.getBIC()
       print thismodel+" results: chisq, K, N, BIC =",LT.minchisq,LT.K,LT.N,BIC[thismodel]
       
       # Write out simple one-line parameter catalog:
       LT.write_catalog(args.outfile)
       print thismodel+" parameter values written to: "+args.outfile

   # -------------------------------------------------------------------------
   
   # # Make some decision about the nature of this system.
   # 
   # if len(models) > 1:
   # # Compare models and report:
   #     print "BIC = ",BIC
   #     print "Hypothesis test result: Bayes factor in favour of nebula is exp[",-0.5*(BIC['Nebula'] - BIC['Lens']),"]"
   #     print "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"

   # -------------------------------------------------------------------------
   
   print "LensTractor stopping."
      
   return  

# ============================================================================

if __name__ == '__main__':
   
   main()

