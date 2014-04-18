'''
This file is part of the lenstractor project.
Copyright 2012 David W. Hogg (NYU) and Phil Marshall (Oxford).

Description
-----------

Wrapper class for tractor operation, to enable high level options (ie sampling or optimization).

To-do
-----
- debug, test

'''

import numpy as np

import tractor
import lenstractor

# import emcee
emcee_defaults = {}

# ============================================================================

class LensTractor():
    '''
    PURPOSE
      Optimize or sample a Lens or Nebula model using the Tractor.

    COMMENTS

    INPUTS
      data             A list of tractor.images
      model            A model (including a list of tractor.sources)
      by               A mode of operation ['optimizing','sampling']
      using            The required settings
      survey           The name of the survey being worked on

    OUTPUTS

    BUGS

    HISTORY
      2014-04-17       Started Marshall & Agnello (UCSB)
    '''
# ----------------------------------------------------------------------------
    
    def __init__(self,dataset,model,survey,vb=0,noplots=True):
    
        self.name = 'LensTractor'
        self.survey = survey
        self.model = model
        self.vb = vb
        self.noplots = noplots
        
        self.chug = tractor.Tractor(dataset)
        for src in self.model.srcs:
           self.chug.addSource(src)

        return None

# ----------------------------------------------------------------------------
    
    def drive(self,by='sampling',using=emcee_defaults):
        
        self.method = by
        self.settings = using
        
        if self.method == 'sampling':
        
            self.sample()
            
        else:
         
            self.optimize()
        
        return None

# ----------------------------------------------------------------------------
    
    def optimize(self):

        # Do a fit by maximizing the posterior PDF ("optimizing")

        Nrounds = self.settings['Nr']
        Nsteps_optimizing_catalog = self.settings['Nc']
        Nsteps_optimizing_PSFs = self.settings['Np']

        if self.vb: 
           print "Optimizing model:"
           print "   - no. of iterations per round to be spent on catalog: ",Nsteps_optimizing_catalog
           print "   - no. of iterations per round to be spent on PSFs: ",Nsteps_optimizing_PSFs
           print "   - no. of rounds: ",Nrounds

        k = 0
        for round in range(Nrounds):

            if self.vb: print "Fitting "+self.model+": seconds out, round",round

            # Freeze the PSF, sky and photocal, leaving the sources:
            if self.vb: print "Thawing catalog..."
            self.chug.thawParam('catalog')
            for image in self.chug.getImages():
                if self.vb: print "Thawing sky..."
                image.thawParams('sky')
                if self.vb: print "Freezing photocal, WCS, PSF..."
                image.freezeParams('photocal', 'wcs', 'psf')
            if self.vb: 
                print "Fitting "+self.model+": Catalog parameters to be optimized are:",self.chug.getParamNames()
                print "Fitting "+self.model+": Initial values are:",self.chug.getParams()
                print "Fitting "+self.model+": Step sizes:",self.chug.getStepSizes()

            # Optimize sources with initial PSF:
            for i in range(Nsteps_optimizing_catalog):
               dlnp,X,a = self.chug.optimize(damp=3,shared_params=False)
               # print "Fitting "+model+": at step",k,"parameter values are:",chug.getParams()
               if self.vb: 
                   print "Progress: k,dlnp = ",k,dlnp
                   print ""
                   print "Catalog parameters:",self.chug.getParamNames()
                   print "Catalog values:",self.chug.getParams()
               if dlnp == 0: 
                   print "Converged? Exiting..."
                   # Although this only leaves *this* loop...
                   break
               k += 1
            if not self.noplots: lenstractor.Plot_state(
                self.chug,
                self.model.name+'_progress_optimizing_step-%02d_catalog'%k,
                SURVEY=self.survey)

            if Nsteps_optimizing_PSFs > 0:
                # Freeze the sources and sky and thaw the psfs:
                if self.vb: print "Freezing catalog..."
                self.chug.freezeParam('catalog')
                for image in self.chug.getImages():
                    if self.vb: print "Thawing PSF..."
                    image.thawParams('psf')
                    if self.vb: print "Freezing photocal, WCS, sky..."
                    image.freezeParams('photocal', 'wcs', 'sky')
                if self.vb: 
                    print "Fitting PSF: After thawing, zeroth PSF = ",self.chug.getImage(0).psf
                    print "Fitting PSF: PSF parameters to be optimized are:",self.chug.getParamNames()
                    print "Fitting PSF: Initial values are:",self.chug.getParams()
                    print "Fitting PSF: Step sizes:",self.chug.getStepSizes()

                # Optimize everything that is not frozen:
                for i in range(Nsteps_optimizing_PSFs):
                   dlnp,X,a = self.chug.optimize(shared_params=False)
                   if self.vb: print "Fitting PSF: at step",k,"parameter values are:",self.chug.getParams()
                   k += 1
                if self.vb: print "Fitting PSF: After optimizing, zeroth PSF = ",self.chug.getImage(0).psf
                if not self.noplots: lenstractor.Plot_state(
                    self.chug,
                    self.model.name+'_progress_optimizing_step-%02d_catalog'%k,
                    SURVEY=self.survey)

                # BUG: PSF not being optimized correctly - missing derivatives?

        # All rounds done! Plot state:
        lenstractor.Plot_state(
            self.chug,
            self.model.name+'_progress_optimizing_zcomplete',
            SURVEY=self.survey)


        return None
        
# ----------------------------------------------------------------------------
    
    def sample(self):

       return None
       
#              if self.vb: 
#                 print "Sampling model parameters with emcee:"
# 
#              # Freeze the wcs and photocal, leaving the PSFs, sky and sources:
#              for image in chug.getImages():
#                 image.freezeParams('photocal', 'wcs')
#                 # Temporary expt:
#                 image.freezeParams('psf')
# 
#              # Get the thawed parameters:
#              p0 = np.array(chug.getParams())
#              print 'Tractor parameters:'
#              for i,parname in enumerate(chug.getParamNames()):
#                    print '  ', parname, '=', p0[i]
#              ndim = len(p0)
#              print 'Number of parameter space dimensions: ',ndim
# 
#              # Make an emcee sampler that uses our tractor to compute its logprob:
#              nw = 8*ndim
#              sampler = emcee.EnsembleSampler(nw, ndim, chug, threads=4)
# 
#              # Start the walkers off near the initialisation point - 
#              # We need it to be ~1 pixel in position, and not too much
#              # flux restrction... 
# 
#              if modeltype=='Lens':
#                 # The following gets us 0.2" in dec:
#                 psteps = np.zeros_like(p0) + 0.00004
#                 # This could be optimized, to allow more initial freedom in eg flux.
# 
#              elif modeltype=='Nebula':
#                 # Good first guess should be some fraction of the optimization step sizes:
#                 psteps = 0.2*np.array(chug.getStepSizes())
# 
#              # BUG - nebula+lens workflow not yet enabled!
# 
# 
#              print "Initial size (in each dimension) of sample ball = ",psteps
# 
#              pp = emcee.EnsembleSampler.sampleBall(p0, psteps, nw)
#              rstate = None
#              lnp = None
# 
#              # Take a few steps - memory leaks fast! (~10Mb per sec)
#              for step in range(10):
# 
#                    print 'EMCEE: Run MCMC step set:', step
#                    t0 = tractor.Time()
#                    pp,lnp,rstate = sampler.run_mcmc(pp, 50, lnprob0=lnp, rstate0=rstate)
# 
#                    print 'EMCEE: Mean acceptance fraction after', sampler.iterations, 'iterations =',np.mean(sampler.acceptance_fraction)
#                    t_mcmc = (tractor.Time() - t0)
#                    print 'EMCEE: Runtime:', t_mcmc
# 
#                    # Find the current posterior means:
#                    # pbar = np.mean(pp,axis=0)
#                    # print "Mean parameters: ",pbar,np.mean(lnp)
# 
#                    # Find the current best sample:
#                    maxlnp = np.max(lnp)
#                    best = np.where(lnp == maxlnp)
#                    pbest = np.ravel(pp[best,:])
#                    print "EMCEE: Best parameters: ",maxlnp,pbest
#                    chug.setParams(pbest)
#                    chisq = -2.0*chug.getLogLikelihood()
#                    print "EMCEE: chisq at Best pt: ",chisq
#                    if not self.noplots: 
#                       chug.setParams(pbest)
#                       lenstractor.Plot_state(
#                           chug,
#                           model+'_progress_sampling_step-%02d'%step,
#                           SURVEY=args.survey)
# 
# 
#              # Take the last best sample and call it a result:
#              chug.setParams(pbest)
# 
#              print 'EMCEE: Best lnprob, chisq:', maxlnp,chisq
#              # print 'dlnprobs:', ', '.join(['%.1f' % d for d in lnp - np.max(lnp)])
#              # print 'MCMC took', t_mcmc, 'sec'
# 
#              # Make the final plot:
#              lenstractor.Plot_state(
#                  chug,
#                  model+'_progress_sampling_zcomplete',
#                  SURVEY=args.survey)
# 
#         return None
# 
# ============================================================================

if __name__ == '__main__':

   pass
