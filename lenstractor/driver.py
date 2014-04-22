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
import os,subprocess

import tractor
import lenstractor

import emcee
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
        
        self.bestpars = None
        self.maxlnp = None
        self.minchisq = None
        
        self.chug = tractor.Tractor(dataset)
        for src in self.model.srcs:
           self.chug.addSource(src)

        # Freeze the PSFs, wcs and photocal, leaving the sky and sources:
        self.chug.thawParam('catalog')
        for image in self.chug.getImages():
           image.thawParams('sky')
           image.freezeParams('photocal')
           image.freezeParams('wcs')
           image.freezeParams('psf')
           
        return None

# ----------------------------------------------------------------------------
    
    def drive(self,by='sampling',using=emcee_defaults):
        
        self.method = by
        self.settings = using
        
        if self.method == 'sampling':
        
            self.sample()
            
        else:
         
            self.optimize()
        
        self.getBIC()
        
        return None

# ----------------------------------------------------------------------------
# Fit the model to the image data by maximizing the posterior PDF ("optimizing")
 
    def optimize(self):

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

            if self.vb: print "Fitting "+self.model.name+": seconds out, round",round

            if self.vb: 
                print "Fitting "+self.model.name+": Catalog parameters to be optimized are:",self.chug.getParamNames()
                print "Fitting "+self.model.name+": Initial values are:",self.chug.getParams()
                print "Fitting "+self.model.name+": Step sizes:",self.chug.getStepSizes()

            # Optimize sources:
            for i in range(Nsteps_optimizing_catalog):
               dlnp,X,a = self.chug.optimize(damp=3,shared_params=False)
               # print "Fitting "+self.model.name+": at step",k,"parameter values are:",self.chug.getParams()
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

#             if Nsteps_optimizing_PSFs > 0:
#                 # Freeze the sources and sky and thaw the psfs:
#                 if self.vb: print "Freezing catalog..."
#                 self.chug.freezeParam('catalog')
#                 for image in self.chug.getImages():
#                     if self.vb: print "Thawing PSF..."
#                     image.thawParams('psf')
#                     if self.vb: print "Freezing photocal, WCS, sky..."
#                     image.freezeParams('photocal', 'wcs', 'sky')
#                 if self.vb: 
#                     print "Fitting PSF: After thawing, zeroth PSF = ",self.chug.getImage(0).psf
#                     print "Fitting PSF: PSF parameters to be optimized are:",self.chug.getParamNames()
#                     print "Fitting PSF: Initial values are:",self.chug.getParams()
#                     print "Fitting PSF: Step sizes:",self.chug.getStepSizes()
# 
#                 # Optimize everything that is not frozen:
#                 for i in range(Nsteps_optimizing_PSFs):
#                    dlnp,X,a = self.chug.optimize(shared_params=False)
#                    if self.vb: print "Fitting PSF: at step",k,"parameter values are:",self.chug.getParams()
#                    k += 1
#                 if self.vb: print "Fitting PSF: After optimizing, zeroth PSF = ",self.chug.getImage(0).psf
#                 if not self.noplots: lenstractor.Plot_state(
#                     self.chug,
#                     self.model.name+'_progress_optimizing_step-%02d_catalog'%k,
#                     SURVEY=self.survey)
# 
#                 # BUG: PSF not being optimized correctly - missing derivatives?

        # Save the best parameters!

        self.maxlnp = self.chug.getLogProb()
        self.bestpars = self.chug.getParams()
        if self.vb: print "Optimizer: Best parameters: ",self.maxlnp,self.bestpars

        self.minchisq = -2.0*self.chug.getLogLikelihood()
        if self.vb: print "Optimizer: chisq at highest lnprob point: ",self.minchisq

        # All rounds done! Plot state:
        lenstractor.Plot_state(
            self.chug,
            self.model.name+'_progress_optimizing_zcomplete',
            SURVEY=self.survey)


        return None
        
# ----------------------------------------------------------------------------
# Compare the model to the image data by sampling the posterior PDF

    def sample(self):

        if self.vb: 
           print "Sampling model parameters with emcee:"

        # Magic numbers:
        nwalkers_per_dim = self.settings['nwp']
        nsnapshots = self.settings['ns']
        nsteps_per_snapshot = self.settings['nss']
        

        # Get the thawed parameters:
        p0 = np.array(self.chug.getParams())
        if self.vb: print 'Tractor parameters:'
        for i,parname in enumerate(self.chug.getParamNames()):
              print '  ', parname, '=', p0[i]
        ndim = len(p0)
        if self.vb: print 'Number of parameter space dimensions: ',ndim

        # Make an emcee sampler that uses our tractor to compute its logprob:
        nw = nwalkers_per_dim*ndim # 8*ndim
        sampler = emcee.EnsembleSampler(nw, ndim, self.chug, threads=4)

        # Start the walkers off near the initialisation point - 
        # We need it to be ~1 pixel in position, and not too much
        # flux restriction... 

        if self.model.name=='Lens':
           # The following gets us 0.2" in dec:
           psteps = np.zeros_like(p0) + 0.00004
           # This could be optimized, to allow more initial freedom in eg flux.

        else:
           # Good first guess should be some fraction of the optimization step sizes:
           psteps = 0.2*np.array(self.chug.getStepSizes())

        # BUG - nebula+lens workflow not yet enabled!
        # psteps should really be a property of the model...


        if self.vb: print "Initial size (in each dimension) of sample ball = ",psteps

#        pp = emcee.EnsembleSampler.sampleBall(p0, psteps, nw)
        pp = emcee.utils.sample_ball(p0, psteps, nw)

        rstate = None
        lnp = None

        # Take a few steps - memory leaks fast! (~10Mb per sec)
        for snapshot in range(nsnapshots):

              if self.vb: print 'Emcee: MCMC snapshot:', snapshot
              t0 = tractor.Time()
              pp,lnp,rstate = sampler.run_mcmc(pp, nsteps_per_snapshot, lnprob0=lnp, rstate0=rstate)

              if self.vb: print 'Emcee: Mean acceptance fraction after', sampler.iterations, 'iterations =',np.mean(sampler.acceptance_fraction)
              t_mcmc = (tractor.Time() - t0)
              if self.vb: print 'Emcee: Runtime:', t_mcmc

              # Find the current posterior means:
              # pbar = np.mean(pp,axis=0)
              # print "Mean parameters: ",pbar,np.mean(lnp)

              # Find the current best sample:
              self.maxlnp = np.max(lnp)
              best = np.where(lnp == self.maxlnp)
              self.bestpars = np.ravel(pp[best,:])
              
              if self.vb: print "Emcee: Best parameters: ",self.maxlnp,self.bestpars
              self.chug.setParams(self.bestpars)
              
              self.minchisq = -2.0*self.chug.getLogLikelihood()
              if self.vb: print "Emcee: chisq at highest lnprob point: ",self.minchisq
              if not self.noplots: 
                 self.chug.setParams(self.bestpars)
                 lenstractor.Plot_state(
                     self.chug,
                     self.model.name+'_progress_sampling_snapshot-%02d'%snapshot,
                     SURVEY=self.survey)


        if self.vb: print 'Emcee: total run time', t_mcmc, 'sec'

        # Make the final plot:
        lenstractor.Plot_state(
            self.chug,
            self.model.name+'_progress_sampling_zcomplete',
            SURVEY=self.survey)

        return None
# ----------------------------------------------------------------------------

    def getBIC(self):    

       self.K = len(self.bestpars)
       self.N = self.chug.getNdata()
       self.BIC = self.minchisq + self.K*np.log(1.0*self.N)
       
       return self.BIC

# ----------------------------------------------------------------------------

    def write_catalog(self,outfile):    

        # Get parameter names and values:
        parnames = self.chug.getParamNames()

        values = np.array(np.outer(1,self.bestpars))

        # Get image names:
        imgnames = []
        for image in self.chug.getImages():
            imgnames.append(image.name)

        # Open up a new file, over-writing any old one:
        try: os.remove(outfile)
        except OSError: pass
        output = open(outfile,'w')

        # Write header:
        hdr = []
        hdr.append('# LensTractor output parameter catalog')
        # hdr.append('# ')
        # hdr.append('# Date: %s' % datestring)
        hdr.append('# ')
        hdr.append('# Model: %s' % self.model.name)
        hdr.append('# Notes:')
        hdr.append('# * First source is always the galaxy, point sources follow')
        for ii,imgname in enumerate(imgnames):
            hdr.append('# * images.image%d = %s' % (ii,imgname))
        hdr.append('# ')
        # Last line contains the parameter names:
        nameline = "#  "
        for name in parnames:
            nameline += name+"  "
        hdr.append(nameline)
        # Write to file:
        for line in hdr:
            output.write("%s\n" % line)    
        # Close file:
        output.close()

        np.savetxt('junk', values)
        cat = subprocess.call("cat junk >> " + outfile, shell=True)
        rm = subprocess.call("rm junk", shell=True)
        if cat != 0 or rm != 0:
          print "Error: write subprocesses failed in some way :-/"
          sys.exit()

        return

           
# ============================================================================

if __name__ == '__main__':

   pass
