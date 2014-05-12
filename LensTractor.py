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
     --manual                  catalog Initialize model positions from catalog
     --optimization-rounds        Nr   Number of rounds of optimization [2]
     --optimization-steps-catalog Nc   Number of steps per round spent
                                        optimizing source catalog [10]
     --optimization-steps-psf     Np   Number of steps per round spent
                                        optimizing PSF catalog [2]
     -o --output             outstem   Stem of output catalog filename       

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
       -o examples/ps1/H1413+117_10x10arcsec \
          examples/ps1/H1413+117_10x10arcsec_55*fits > \
          examples/ps1/H1413+117_10x10arcsec_Nebula4.log
       
     python LensTractor.py -n 2 \
       -o examples/sdss/0951+2635/0951+2635 \
          examples/sdss/0951+2635/*sci.fits > \
          examples/sdss/0951+2635/0951+2635_Nebula2.log

     set id = KIDS_SLID_10058881_SID_8668
     python LensTractor.py -v -l -z --survey KIDS \
       -o examples/kids/${id} \
          examples/kids/${id}_u_???.fits \
          examples/kids/${id}_g_???.fits \
          examples/kids/${id}_r_???.fits > \
          examples/kids/${id}.log &

   
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
   # Optimizing only:
   parser.add_argument('-z', '--optimize', dest='optimize', action='store_true', default=False, help='Optimize posterior PDF')
   # Sampling only:
   parser.add_argument('-s', '--sample', dest='MCMC', action='store_true', default=False, help='Sample posterior PDF')
   # Plotting:
   parser.add_argument('-x', '--no-plots', dest='noplots', action='store_true', default=False, help='Skip plotting')
   # Lens model only:
   parser.add_argument('-l', '--lens', dest='lens', action='store_true', default=False, help='Fit Lens model')
   # Nebula model only:
   parser.add_argument('-n', '--nebula', dest='K', type=int, default=0, help='Fit NebulaK model, provide K')
   # Output filename:
   parser.add_argument('-o', '--output', dest='outstem', type=str, default='lenstractor.cat', help='Output catalog filename stem')
   # Survey we are working on (affects data read-in):
   parser.add_argument('--survey', dest='survey', type=str, default="PS1", help="Survey (SDSS, PS1 or KIDS)")
   # Manual input of model initialization:
   parser.add_argument('--manual', dest='catalog', type=str, default="None", help="Catalog of Nebula model parameters, for initializing positions")

   # Read in options and arguments - note only sci and wht images are supplied:
   args = parser.parse_args()
   
      
   if (len(args.inputfiles) < 1):
      # parser.print_help()
      print main.__doc__  # Whoah! What does this do?! Some sort of magic.
      sys.exit(-1)
   
   vb = args.verbose
   
   # Manual initialization:
   if args.catalog != 'None':
      manual = True
   else:
      manual = False
      
   # Workflow:
   if args.lens:
      # modelnames = ['Nebula2','Lens']
      modelnames = ['Nebula4','Lens']
   elif args.K > 0:
      modelnames = ['Nebula'+str(args.K)]
   else:
      modelnames = ['Nebula2','Nebula4','Lens']
         
   # BIC = dict(zip(modelnames,np.zeros(len(modelnames))))
   BIC = dict()


   if vb: 
      print "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
      print "                               LensTractor "
      print "    Fitting",modelnames," models to a deck of FITS postcards"
      print "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
 
   # -------------------------------------------------------------------------
   # Read in images (using IO functions in dm.py)
   
   # Organise the deck of inputfiles into scifiles and varfiles:
   scifiles,varfiles = lenstractor.Riffle(args.inputfiles,vb=vb)
   
   # Read into Tractor Image objects, and see what filters we have:   
   images,centroids,total_mags,bands = lenstractor.Deal(scifiles,varfiles,SURVEY=args.survey,vb=vb)
   
   # Estimate object centroid and SED:
   position,SED = lenstractor.Turnover(bands,total_mags,centroids,vb=vb)
   
   # Package up:
   dataset = list(images)
            
   # -------------------------------------------------------------------------
   # Step through all the models in the workflow, initialising and fitting:
   
   previous = None
   counter = 0
   for modelname in modelnames: 
      
       if vb: 
           print "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
           print "Initializing model: "+modelname       
                     
       model = lenstractor.Model(modelname,vb=vb)
       
       if modelname != 'Lens':
          if manual:
              model.initialize('from_scratch', position=args.catalog, SED=SED)
          
          else:
              if previous is None:
                  model.initialize('from_scratch', position=position, SED=SED)
              else:
                  model.initialize(previous)
       else:
          model.initialize(previous)
       
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
       
       LT = lenstractor.LensTractor(dataset,model,args.outstem,args.survey,counter=counter,vb=vb,noplots=args.noplots)

       # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

       if args.optimize:

           LT.drive(by='optimizing')

       elif args.MCMC:

           LT.drive(by='sampling')
           
       else:

           LT.drive(by='cunning_and_guile')
           
       # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       
       if vb:
           print "Fit complete."
           print "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
       
       # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              
       # Compute BIC for this fit:
       BIC[model.flavor] = LT.getBIC()
       if vb: print modelname+" results: chisq, K, N, BIC =",LT.minchisq,LT.K,LT.N,BIC[model.flavor]
       
       # Have model print itself:
       for component in model.srcs:
           if vb: print component
       
       # Write out simple one-line parameter catalog:
       outfile = LT.write_catalog()
       if vb: print modelname+" parameter values written to: "+outfile

       # Save Nebula2 or Nebula4? Depends on BIC...
       previous = model.copy()
       counter = LT.counter + 1

   # -------------------------------------------------------------------------
   
   # Make some decision about the nature of this system.
   
   if len(modelnames) > 1:
       assert model.name == 'Lens'
       
   # Compare models and report:
       if vb: print "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
       if vb: print "BIC = ",BIC
       logBF = -0.5*(BIC['Lens'] - BIC['Nebula'])
       Nimg = model.srcs[0].getMultiplicity()
       if vb: print "Hypothesis test result: Bayes factor in favour of Lens is exp[",logBF,"]"
       if Nimg < 2:
           if vb: print "BUT: Lens predicts only 1 image, so it's not a strong lens."
       if vb: print "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"

       # Record result as a "cookie":
       if logBF > 0 and Nimg > 1:
           result = 'Lens'
       else:
           result = 'Nebula'
       cookie = LT.set_cookie(args.outstem, result)
       if vb: print 'Result "'+result+'" written to '+cookie
       
   # -------------------------------------------------------------------------
   
   print "LensTractor stopping."
      
   return  

# ============================================================================

if __name__ == '__main__':
   
   main()

