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

import pylab as plt

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
      plot             Make plots if desired

    OUTPUTS

    BUGS

    HISTORY
      2014-04-17       Started Marshall & Agnello (UCSB)
    '''
# ----------------------------------------------------------------------------
    
    def __init__(self,dataset,model,outstem,survey,counter=0,vb=0,noplots=True):
    
        self.name = 'LensTractor'
        self.survey = survey

        self.settings = {}
        # Optimization settings:
 #       self.settings['Nrounds'] = 5
 #       self.settings['Nsteps_optimizing_catalog'] = 100
 #       self.settings['Nsteps_optimizing_PSFs'] = 2
        # Sampling settings:
        self.settings['Nwalkers_per_dim'] = 8
        self.settings['Nsnapshots'] = 3
        self.settings['Nsteps_per_snapshot'] = 5000
        self.settings['Restart'] = True
        
        self.model = model
        
        self.vb = vb
        self.noplots = noplots
        self.plot_all = True
        
        self.outstem = outstem.split('.')[0]
        
        self.bestpars = None
        self.maxlnp = None
        self.minchisq = None
        self.psteps = None
        
        self.counter = counter
        
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
        
        # Plot initial state:
        if not self.noplots: 
            self.plot_state('progress-%02d_initial_'%self.counter+self.model.name)
        
        return None

# ----------------------------------------------------------------------------
# Drive the LensTractor. We have both steepest ascent and MCMC capability.
# Try a mixture!
   
    def drive(self,by='cunning_and_guile'):
        
        self.method = by
        
        if self.method == 'sampling':
        
            self.sample()
            
        elif self.method == 'optimizing':
         
            # First optimize to get the model about right, at fixed PSF:
            self.settings['Nrounds'] = 1
            self.settings['Nsteps_optimizing_catalog'] = 100000
            self.settings['Nsteps_optimizing_PSFs'] = 0
            self.optimize()
        
            # Now optimize PSF at fixed model:
            self.settings['Nrounds'] = 1
            self.settings['Nsteps_optimizing_catalog'] = 0
            self.settings['Nsteps_optimizing_PSFs'] = 2
            self.optimize()
            
            # Refine model at best PSF:
            self.settings['Nrounds'] = 1
            self.settings['Nsteps_optimizing_catalog'] = 10000
            self.settings['Nsteps_optimizing_PSFs'] = 0
            self.optimize()
        
        
        else: # Apply cunning and guile!
         
            # First optimize to get the fluxes about right:
            self.settings['Nrounds'] = 1
            self.settings['Nsteps_optimizing_catalog'] = 100000
            self.settings['Nsteps_optimizing_PSFs'] = 0
            self.optimize()          
            
            # Now draw a few samples to shuffle the positions:
            self.settings['Nsnapshots'] = 3
            self.settings['Nwalkers_per_dim'] = 2
            self.settings['Nsteps_per_snapshot'] = 4000
            self.settings['Restart'] = True
            self.sample()
            
            # Now optimize to refine model and PSF:
            self.settings['Nrounds'] = 1
            self.settings['Nsteps_optimizing_catalog'] = 50000
            self.settings['Nsteps_optimizing_PSFs'] = 0
            self.optimize()
            self.settings['Nrounds'] = 1
            self.settings['Nsteps_optimizing_catalog'] = 0
            self.settings['Nsteps_optimizing_PSFs'] = 2
            self.optimize()
            self.settings['Nrounds'] = 1
            self.settings['Nsteps_optimizing_catalog'] = 10000
            self.settings['Nsteps_optimizing_PSFs'] = 0
            self.optimize()
            
        self.getBIC()
        
        return None

# ----------------------------------------------------------------------------
# Fit the model to the image data by maximizing the posterior PDF 
# ("optimizing") with respect to the parameters.
  
    def optimize(self):

        Nrounds = self.settings['Nrounds']
        Nsteps_optimizing_catalog = self.settings['Nsteps_optimizing_catalog']
        Nsteps_optimizing_PSFs = self.settings['Nsteps_optimizing_PSFs']

        if self.vb: 
           print " "
           print "Optimizing model:"
           print "   - no. of iterations per round to be spent on catalog: ",Nsteps_optimizing_catalog
           print "   - no. of iterations per round to be spent on PSFs: ",Nsteps_optimizing_PSFs
           print "   - no. of rounds: ",Nrounds

        for round in range(Nrounds):
            self.counter += 1

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
                   print "Progress: counter,dlnp = ",self.counter,dlnp
                   print ""
                   print "Catalog parameters:",self.chug.getParamNames()
                   print "Catalog values:",self.chug.getParams()
               if dlnp == 0: 
                   print "Converged? Exiting..."
                   # Although this only leaves *this* loop...
                   break
            
            if not self.noplots: 
                self.plot_state('progress-%02d_optimizing_'%self.counter+self.model.name)

            if Nsteps_optimizing_PSFs > 0:
            
                # Freeze the sources and calibration, and thaw the psfs:
                if self.vb: print "Freezing catalog..."
                self.chug.freezeParams('catalog')
                for image in self.chug.getImages():
                    if self.vb: print "Thawing PSF..."
                    image.thawParams('psf')
                    if self.vb: print "Freezing photocal, WCS, sky (just to make sure...)"
                    image.freezeParams('photocal', 'wcs', 'sky')
                if self.vb: 
                    print "Fitting PSF: After thawing, zeroth PSF = ",self.chug.getImage(0).psf
                    print "Fitting PSF: PSF parameters to be optimized are:",self.chug.getParamNames()
                    print "Fitting PSF: Initial values are:",self.chug.getParams()
                    print "Fitting PSF: Step sizes:",self.chug.getStepSizes()

                # Optimize everything that is not frozen:
                for i in range(Nsteps_optimizing_PSFs):
                   dlnp,X,a = self.chug.optimize(shared_params=False)
                   if self.vb: print "Fitting PSF: at counter =",self.counter,"parameter values are:",self.chug.getParams()
                self.counter += 1
                
                if self.vb: print "Fitting PSF: After optimizing, zeroth PSF = ",self.chug.getImage(0).psf
                if not self.noplots: self.plot_state(
                    'progress-%02d_optimizing_PSF_for_'%self.counter+self.model.name)

                # Freeze the psf again, and thaw the sources:
                if self.vb: print "Re-thawing catalog..."
                self.chug.thawParams('catalog')
                for image in self.chug.getImages():
                    if self.vb: print "Re-freezing PSF..."
                    image.freezeParams('psf')
                    if self.vb: print "Re-freezing photocal, WCS, sky (just to make sure...)"
                    image.freezeParams('photocal', 'wcs', 'sky')


        # Save the best parameters!

        self.maxlnp = self.chug.getLogProb()
        self.bestpars = self.chug.getParams()
        if self.vb: print "Optimizer: Best parameters: ",self.maxlnp,self.bestpars

        self.minchisq = -2.0*self.chug.getLogLikelihood()
        if self.vb: print "Optimizer: chisq at highest lnprob point: ",self.minchisq

        return None
        
# ----------------------------------------------------------------------------
# Fit the model to the image data by drawing samples from the posterior PDF 
# that have high probability density: note, this is not really sampling,
# its *sampling to optimize*...

    def sample(self):

        if self.vb: 
           print "Sampling model parameters with Emcee:"

        # Magic numbers:
        Nwalkers_per_dim = self.settings['Nwalkers_per_dim']
        Nsnapshots = self.settings['Nsnapshots']
        Nsteps_per_snapshot = self.settings['Nsteps_per_snapshot']
        Restart = self.settings['Restart']

        # Get the thawed parameters:
        p0 = np.array(self.chug.getParams())
        if self.vb: print 'Tractor parameters:'
        for i,parname in enumerate(self.chug.getParamNames()):
              print '  ', parname, '=', p0[i]
        Ndim = len(p0)
        if self.vb: print 'Number of parameter space dimensions: ',Ndim

        # Make an emcee sampler that uses our tractor to compute its logprob:
        Nw = Nwalkers_per_dim*Ndim # 8*ndim
        sampler = emcee.EnsembleSampler(Nw, Ndim, self.chug, threads=4)

        # Start the walkers off near the initialisation point - 
        # We need it to be ~1 pixel in position, and not too much
        # flux restriction... But use any psteps we already have!
        
        if self.psteps is None:
            if self.model.name=='Lens':
               # The following gets us 0.05" in dec:
               self.psteps = np.zeros_like(p0) + 0.00001
               # This could be optimized, to allow more initial freedom in eg flux.
            else:
               # Good first guess should be some fraction of the optimization step sizes:
               self.psteps = 0.01*np.array(self.chug.getStepSizes())
               
        if self.vb: print "Initial size (in each dimension) of sample ball = ",self.psteps
        
        #pp = emcee.EnsembleSampler.sampleBall(p0, self.psteps, Nw)
        pp = emcee.utils.sample_ball(p0, self.psteps, Nw)
        rstate = None
        lnp = None

        # Take a few steps - memory leaks fast! (~10Mb per sec)
        for snapshot in range(1,Nsnapshots+1):
              self.counter += 1
              
              if self.vb: print 'Emcee: MCMC snapshot:', snapshot
              t0 = tractor.Time()
              pp,lnp,rstate = sampler.run_mcmc(pp, Nsteps_per_snapshot, lnprob0=lnp, rstate0=rstate)

              if self.vb: print 'Emcee: Mean acceptance fraction after', sampler.iterations, 'iterations =',np.mean(sampler.acceptance_fraction)
              t_mcmc = (tractor.Time() - t0)
              if self.vb: print 'Emcee: Runtime:', t_mcmc

              # Find the current best sample, and sample ball:
              self.maxlnp = np.max(lnp)
              best = np.where(lnp == self.maxlnp)
              self.bestpars = np.ravel(pp[best,:])
              if self.vb: print "Emcee: Best parameters: ",self.maxlnp,self.bestpars
              
              self.minchisq = -2.0*self.chug.getLogLikelihood()
              if self.vb: print "Emcee: chisq at highest lnprob point: ",self.minchisq
              if not self.noplots: 
                 self.chug.setParams(self.bestpars)
                 self.plot_state('progress-%02d_sampling_'%self.counter+self.model.name)

              if Restart:
                 # Make a new sample ball centred on the current best point,
                 # and with width given by the standard deviations in each
                 # dimension:
                 self.chug.setParams(self.bestpars)
                 p0 = np.array(self.chug.getParams())
                 self.psteps = np.std(pp,axis=0)
                 # pp = emcee.EnsembleSampler.sampleBall(p0, self.psteps, Nw)
                 pp = emcee.utils.sample_ball(p0, self.psteps, Nw)
                 rstate = None
                 lnp = None


        if self.vb: print 'Emcee: total run time', t_mcmc, 'sec'


        return None
        
# ----------------------------------------------------------------------------

    def getBIC(self):    

       self.K = len(self.bestpars)
       self.N = self.chug.getNdata()
       self.BIC = self.minchisq + self.K*np.log(1.0*self.N)
       
       return self.BIC

# ----------------------------------------------------------------------------

    def write_catalog(self):    

        # Get parameter names and values:
        parnames = self.chug.getParamNames()

        values = np.array(np.outer(1,self.bestpars))

        # Get image names:
        imgnames = []
        for image in self.chug.getImages():
            imgnames.append(image.name)

        # Set catalog file name:
        outfile = self.outstem+'_'+self.model.name+'.cat'
        
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

        return outfile

           
# ----------------------------------------------------------------------------

    def set_cookie(self,outstem,result):    

        outfile = self.outstem+'_result.cookie'
        
        # Open up a new file, over-writing any old one:
        try: os.remove(outfile)
        except OSError: pass
        output = open(outfile,'w')

        # Write result to file:
        output.write("%s\n" % result)    
        # Close file:
        output.close()

        return outfile

           
# ----------------------------------------------------------------------------
# Plot progress.

    def plot_state(self,suffix):
       '''
       Make all the plots we need to assess the state of the LensTractor. 
       Mainly, a multi-panel figure of image, synthetic image and chi, for 
       each image being modelled.

       self.chug is a Tractor object, containing a list of images.
       '''

       if self.plot_all:
           # Make one giant plot with all progress panels on it, and 
           # save it to PNG:

           # Plotting setup:
           px,py = 4,len(self.chug.images)
           figprops = dict(figsize=(5*px,5*py), dpi=128)
           adjustprops = dict(\
               left=0.05,\
               bottom=0.05,\
               right=0.95,\
               top=0.95,\
               wspace=0.1,\
               hspace=0.1)
           
           # Start the plot:
           fig = plt.figure(**figprops)
           fig.subplots_adjust(**adjustprops)
           plt.clf()
           plt.gray()
           counter = 0
          
           # Name the plot file:
           pngfile = self.outstem+'_'+suffix+'.png'


       else:
           # Make one plot per image, and save each to PNG.

           # Plotting setup:
           px,py = 2,2
           figprops = dict(figsize=(5*px,5*py), dpi=128)
           adjustprops = dict(\
              left=0.05,\
              bottom=0.05,\
              right=0.95,\
              top=0.95,\
              wspace=0.1,\
              hspace=0.1)



       # Loop over images, making plots:

       for i,image in enumerate(self.chug.images):

          if image.name is None:
             imname = suffix+str(i)
          else:
             imname = image.name

          chi = self.chug.getChiImage(i)

          if self.survey == 'PS1':
              ima, chia, psfa = lenstractor.PS1_imshow_settings(image,chi)
          elif self.survey == 'KIDS':
              ima, chia, psfa = lenstractor.KIDS_imshow_settings(image,chi)
          else:
              # Do the same as for PS1
              scale = np.sqrt(np.median(1.0/image.invvar[image.invvar > 0.0]))
              ima = dict(interpolation='nearest', origin='lower',
                         vmin=-100.*scale, vmax=3.*scale)

              chia = dict(interpolation='nearest', origin='lower',
                          vmin=-5., vmax=5.)

              psfa = dict(interpolation='nearest', origin='lower')

          if not self.plot_all:
              # Start a new plot:
              fig = plt.figure(**figprops)
              fig.subplots_adjust(**adjustprops)
              plt.clf()
              plt.gray()
              counter = 0

          # 1) Data Image
          counter += 1
          plt.subplot(py,px,counter)
          plt.imshow(-image.data, **ima)
          self.tidyup_plot()
          plt.title('Observed image')
          # Overlay image filter in lower left corner
          plt.text(1,1,image.photocal.bandname)
          # Figure out how to get this in bottom right hand corner instead

          # 2) Predicted image
          counter += 1
          plt.subplot(py,px,counter)
          model = self.chug.getModelImages()[i]
          plt.imshow(-model, **ima)
          self.tidyup_plot()
          plt.title('Predicted image')
          # Overlay name of model in lower left corner
          plt.text(1,1,self.model.name)
          # Figure out how to get this in top left hand corner instead
          # Overlay cartoon of model...
          # TBI!

          # 3) Normalised residual
          counter += 1
          plt.subplot(py,px,counter)
          plt.imshow(-chi, **chia)
          self.tidyup_plot()
          if self.survey == 'KIDS':
              # It is not clear why the residual image is not in units of
              # sigma. Perhaps this causes problems in the modelling.
              # This code is not refactored into kids.py since it should
              # not be necessary in the first place.
              plt.title('Residuals (flexible scale)')
          else:
              plt.title('Residuals ($\pm 5\sigma$)')
          # Overlay quantified goodness of fit, in sigma from acceptable...
          # TBI!

          # 4) PSF image
          counter += 1
          plt.subplot(py,px,counter)
          psfimage = image.psf.getPointSourcePatch(*model.shape).patch
          plt.imshow(-psfimage, **psfa)
          self.tidyup_plot()
          plt.title('PSF')

          if not self.plot_all:
              # Save this image's plot:
              pngfile = imname+'_'+suffix+'.png'
              plt.savefig(pngfile)
              if self.vb:
                  print "Progress snapshot saved to "+pngfile

       if self.plot_all:
           # Save the giant plot:
           plt.savefig(pngfile)
           if self.vb:
               print "Progress snapshot saved to "+pngfile
       

       return

    # ----------------------------------------------------------------------------
    # Turn off the axis ticks and labels:

    def tidyup_plot(self):
       ax = plt.gca()
       ax.xaxis.set_ticks([])
       ax.yaxis.set_ticks([])
       return

# ============================================================================

if __name__ == '__main__':

   pass
