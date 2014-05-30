'''
This file is part of the LensFinder project.
Copyright 2012 David W. Hogg (NYU) and Phil Marshall (Oxford).

Description
-----------

Simple lens models for lensed quasar detection. Given deflector parameters (position, 
Einstein radius, gammax, phix), and a source position, solve for the 2,3 or 4 images, 
and return image positions and fluxes (and many other things).


To-do
-----
- test lens-equation solving [Done: 99.94% success, 10ms per solve]
- perhaps write unit tests for potentials(), deflections(), magnifications()
- write down priors on magnification / mag bias
'''

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import rc
    rc('font',**{'family':'serif','serif':'Computer Modern Roman','size':18})
    rc('text', usetex=True)

import numpy as np
# import pylab as plt
import matplotlib.pyplot as plt

# Caustic drawing:
from matplotlib import nxutils
#from matplotlib.path import Path  
# Requires matplotlib v1.3.1 or greater; upgrade with 'pip install --upgrade matplotlib'

# Global verbosity:
vb = 0

# ============================================================================

# cheap and nasty singular isothermal sphere model
# initialize with
# - position        : lens position (same units as einsteinradius)
# - einsteinradius  : Einstein radius (some angular units)
# - gammacos2phi    : external shear amplitude times cos(2 phi)
# - gammasin2phi    : external shear amplitude times sin(2 phi)

class GravitationalLens:

    def __init__(self, position, einsteinradius, gamma, phi):
        self.name = 'Gravitational lens, modeled as SIS density profile + external shear'
        self.position = np.array(position)
        self.einsteinradius = einsteinradius
        self.set_gamma_phi(gamma, phi)
        return None

    def __str__(self):
        return '%s (%s %.2f %.2f %.2f)' % (self.name, self.position, self.einsteinradius, self.gamma, self.phi)

# ----------------------------------------------------------------------------
# External shear set-up:

    def set_gamma_trig(self, gammacos2phi, gammasin2phi):
        # self.phi is in radians
        self.gammacos2phi = gammacos2phi
        self.gammasin2phi = gammasin2phi
        self.gamma = np.sqrt(self.gammacos2phi**2 + self.gammasin2phi**2)
        self.phi = 0.5 * np.arctan2(self.gammasin2phi, self.gammacos2phi)
        return None

    def set_gamma_phi(self, gamma, phi):
        # phi is in radians
        self.gamma = gamma
        self.phi = phi
        self.gammacos2phi = self.gamma * np.cos(2. * self.phi)
        self.gammasin2phi = self.gamma * np.sin(2. * self.phi)
        return None

# ----------------------------------------------------------------------------
# Caustics and critical curves
#   magic number 1024 sets resolution (plots are made with points not lines)

    def eigs(self):
   	  return (np.array([np.cos(self.phi), np.sin(self.phi)]), 
                np.array([-np.sin(self.phi), np.cos(self.phi)]))
    
    def cross(self):
        # The sourceposition = 0 morphology; undefined when gamma=0
        assert(self.gamma > 0.)
        eig1, eig2 = self.eigs()
        foo = np.array([ eig1 / (1. - self.gamma),
                         eig2 / (1. + self.gamma),
                         -eig1 / (1. - self.gamma),
                         -eig2 / (1. + self.gamma)])
        return self.position + self.einsteinradius * foo
    
    def critical_curve(self, npts=1024):
        eig1, eig2 = self.eigs()
        ddphi = 2. * np.pi / npts
        dphi = np.arange(0., 2. * np.pi + 0.5 * ddphi, ddphi)
        cosphi = np.cos(dphi)
        sinphi = np.sin(dphi)
        costwophi = np.cos(2.0*dphi)
        r = self.einsteinradius*(1 - self.gamma*costwophi)/(1.0-self.gamma*self.gamma)
        critcurve = np.outer(np.ones_like(dphi), self.position)
        critcurve += np.outer(r*cosphi, eig1)
        critcurve += np.outer(r*sinphi, eig2)
        return critcurve

    def radial_caustic(self, npts=1024):
        # note crazy over-wrapping on dphi; this is used in guess_image_positions()
        # note computational overkill; it's a freakin' circle.
        eig1, eig2 = self.eigs()
        ddphi = 2. * np.pi / npts
        dphi = np.arange(0., 2. * np.pi + 1.5 * ddphi, ddphi)
        cosphi = np.cos(dphi)
        sinphi = np.sin(dphi)
        caustic = np.outer(np.ones_like(dphi), self.position)
        caustic += np.outer(self.einsteinradius * cosphi, eig1)
        caustic += np.outer(self.einsteinradius * sinphi, eig2)
        
        return caustic

    def tangential_caustic(self, npts=1024):
        eig1, eig2 = self.eigs()
        ddphi = 2. * np.pi / npts
        dphi = np.arange(0., 2. * np.pi + 0.5 * ddphi, ddphi)
        sincubedphi = np.sin(dphi)**3
        coscubedphi = np.cos(dphi)**3
        betaplus = 2.0*self.einsteinradius*self.gamma/(1.0+self.gamma)
        betaminus = 2.0*self.einsteinradius*self.gamma/(1.0-self.gamma)
        caustic = np.outer(np.ones_like(dphi), self.position)
        caustic += np.outer(-betaplus*coscubedphi, eig1)
        caustic += np.outer(betaminus*sincubedphi, eig2)
        
        return caustic

# ----------------------------------------------------------------------------
# Expected number of images
#   returns 1, 2, 3, or 4
    
    def number_of_images(self, sourceposition):
        rc=self.radial_caustic()
        tc=self.tangential_caustic()

         # Old Matplotlib:
        if nxutils.points_inside_poly(np.atleast_2d(sourceposition),rc)[0]:
            if nxutils.points_inside_poly(np.atleast_2d(sourceposition),tc)[0]:
                return 4
            return 2
        if nxutils.points_inside_poly(np.atleast_2d(sourceposition),tc)[0]:
            return 3
        
        # New Matplotlib:
        # rc=Path(rc)
        # tc=Path(tc)
        # if rc.contains_points(np.atleast_2d(sourceposition))[0]:
        #    if tc.contains_points(np.atleast_2d(sourceposition))[0]:
        #        return 4
        #    return 2
        # if tc.contains_points(np.atleast_2d(sourceposition))[0]:
        #    return 3

        return 1
        

# ----------------------------------------------------------------------------
# Getting started - initial guesses for image positions
#
#   inputs: a single source position, shape (2)
#   outputs: a first guess at N image positions, shape (N, 2)
#   notes:
#   - N=3 is fugly
#   - N=4 is worse

    def guess_image_positions(self, sourceposition):
        N = self.number_of_images(sourceposition)
        assert(N > 0 and N < 5)
        if N == 1:
            ipos = 2.0*sourceposition.copy()
        if N == 2:
            dpos = sourceposition.copy() - self.position
            # radial vector, length = thetaE
            blah = self.einsteinradius * dpos / np.sqrt(np.dot(dpos, dpos))
            ipos = np.array([sourceposition.copy() + blah, sourceposition.copy() - blah])
        if N == 3:
            dpos = sourceposition.copy() - self.position
            # tangential vector, length = 5*thetaE to cope with very 
            # wide separation arcs...
            blah = 5.0*self.einsteinradius * dpos / np.sqrt(np.dot(dpos, dpos))
            rot = np.array([[0,1],[-1,0]])
            wffl = np.dot(rot,blah)
            ipos = np.array([sourceposition.copy() - wffl, 
                             sourceposition.copy(), 
                             sourceposition.copy() + wffl])
        if N == 4:
            ipos = self.guess_ring_positions(sourceposition)
            ipos = np.array([self.guess_radial_position(sourceposition, ip) for ip in ipos])
        return ipos

    def guess_ring_positions(self, sourceposition):
        # BUG: does not necessarily return 4 images! Hack to cross if fail.
        thisvb = 0
        Nmin, Nmax = 0, 0
        npts = 1024
        while (Nmin != 2 or Nmin != Nmax) and npts < 1e6:
            if thisvb: print 'making ring of',npts
            ringpos = self.radial_caustic(npts=npts)
            ringpos = self.position + (1. + self.gamma) * (ringpos - self.position)
            td = self.time_delays(sourceposition, ringpos)
            minI = ((td[1:-1] < td[0:-2]) *
                    (td[1:-1] < td[2:]))
            Nmin = np.sum(minI)
            maxI = ((td[1:-1] > td[0:-2]) *
                    (td[1:-1] > td[2:]))
            Nmax = np.sum(maxI)
            if thisvb: print 'Nmin, Nmax', Nmin, Nmax
            npts = npts * 4
            # note horrifying offset: exercise to the reader
        ringpos = ringpos[1:-1]
        ipos = np.append(ringpos[minI], ringpos[maxI], axis=0)
        
        if len(ipos) != 4:
          ipos = np.zeros([4,2])
          for i in range(4):
            ipos[i,0] = self.einsteinradius * np.cos(i*0.5*np.pi)
            ipos[i,1] = self.einsteinradius * np.sin(i*0.5*np.pi)
        
        return ipos

    def guess_radial_position(self, sourceposition, imageposition):
        # NB. takes just one image position as input, returns update on output
        rfactor = np.arange(0.5, 2.0, 0.005)
        ipos = np.outer(rfactor, imageposition - self.position) + self.position
        return ipos[np.argmin(self.time_delays(sourceposition, ipos))]

# ----------------------------------------------------------------------------
# Refining image positions during solve process.
#
#   input: a single source position shape (2) and a first guess at N image positions shape (N, 2) 
#   output: N image positions shape (N, 2)                                                        
#   magic numbers: tol, too_many (needed to escape infinite oscillation)                          

    def refine_image_positions(self, sourceposition, guessedimagepositions, tol=1.e-20, alpha=1.0, too_many=1e3):
        ipos = np.atleast_2d(guessedimagepositions)
        parities = self.parities(ipos)
        N = len(ipos[:,0])
        dspos = np.outer(np.ones(N), sourceposition) - self.source_positions(ipos)
        i = 0
        while np.sum(dspos**2) > tol and i < too_many:
            i += 1
            dipos = np.array([np.dot(tens, dsp) for tens, dsp in zip(self.magnification_tensors(ipos), dspos)])
            ipos = ipos + alpha * dipos
            dspos = np.outer(np.ones(N), sourceposition) - self.source_positions(ipos)        
        fail = (np.sum(dspos**2) > tol)
        return ipos,fail

# ----------------------------------------------------------------------------
# When stuck, images need perturbing before refinement proceeds.
#
#   input: a set of (failed) N image positions shape (N, 2)
#   output: N perturbed image positions shape (N, 2)

    def perturb_image_positions(self, imagepositions,sigma):
        ipos = np.atleast_2d(imagepositions)
        N = len(ipos)
        ipos[:,0] = ipos[:,0] + np.random.uniform(-sigma,sigma,size=N)
        ipos[:,1] = ipos[:,1] + np.random.uniform(-sigma,sigma,size=N)
        return ipos

# One sign of being stuck is to return two copies of the same image:
    def duplicate_locations(self,xx,invert=False):
        # input: array of 2, 3 or 4 numbers
        # output: indices of the ones that are the same
        uniq = np.unique(xx)
        if not invert:
          index = np.zeros(len(xx) - len(uniq),dtype=np.int)
          j = 0
          for x in uniq:
            i = np.where(xx == x)
            if len(i[0]) > 1:
              index[j] = i[0][0]
              j += 1        
        else:
          index = np.zeros(len(uniq),dtype=np.int)
          j = 0
          for x in uniq:
            i = np.where(xx == x)
            if len(i[0]) == 1:
              index[j] = i[0][0]
              j += 1
        return index

# ----------------------------------------------------------------------------
# Solve for image positions given source position.
#
#   input: source position, shape (2,)
#   output: image positions, shape (N,2)

    def image_positions(self, sourceposition, keepgoing=True):

        assert(sourceposition.shape == (2, ))
        ipos = self.guess_image_positions(sourceposition)
        N = self.number_of_images(sourceposition)
        initialparities = self.parities(ipos)
                
        done = False
        fail = False
        humph = 1
        while not done:
        
          # Reset starting position every 10th round, perturbing farther each time:
          if keepgoing and (humph % 10) == 0: 
            ipos = self.guess_image_positions(sourceposition)
            ipos = self.perturb_image_positions(ipos,0.01*humph*self.einsteinradius)
            initialparities = self.parities(ipos)
          
          # If failed last round, perturb the images and try again:
          if fail and keepgoing:
            ipos = self.perturb_image_positions(ipos,0.1*self.einsteinradius)
          
          # OK, inch forward:
          fail = False
          ipos,fail = self.refine_image_positions(sourceposition, ipos)

          magnifications = self.magnifications(ipos)
          parities = np.sign(magnifications)

          # Check image parities:
          if (N == 2) or (N == 4):
              if (np.sum(parities) != 0):
                  if vb: 
                    print 'image_positions: WARNING: parities wrong, some images have either merged or collapsed'
                    # print 'image positions:',ipos
                  fail = True
          elif (N == 3):
              if (np.sum(parities - initialparities) != 0):
                  if vb: 
                    print 'image_positions: WARNING: parities wrong, some images have merged or collapsed'
                    # print 'image positions:',ipos
                  fail = True
          
          # Check we have the right number of images:
          if len(ipos[:,0]) > len(set(ipos[:,0])):
              if vb: 
                print 'image_positions: WARNING: x positions wrong, some images have collapsed'
                print 'image positions:',ipos
              fail = True
          if len(ipos[:,1]) > len(set(ipos[:,1])):
              if vb: 
                print 'image_positions: WARNING: y positions wrong, some images have collapsed'
                print 'image positions:',ipos
              fail = True
          
          # Check that magnifications are unique:
          if len(magnifications) > len(set(magnifications)):
              if vb: 
                print 'image_positions: WARNING: magnifications wrong, some images have collapsed'
                print 'image magnifications:',magnifications
              fail = True
              # Find duplicated images and re-position them:
              duplicated = False
              if N==4 and len(set(magnifications))==3:
                duplicated = True
                index = self.duplicate_locations(magnifications)
                ipos[index,0] = 0.0
                ipos[index,1] = 0.0
              elif N==4 and len(set(magnifications))==2:
                duplicated = True
                index = self.duplicate_locations(magnifications)
                nindex = self.duplicate_locations(magnifications,invert=True)
                ipos[index[0],0] = 0.0
                ipos[index[0],1] = 0.0
                ipos[index[1],0] = 0.5*(ipos[nindex[0],0]+ipos[nindex[0],0])
                ipos[index[1],1] = 0.5*(ipos[nindex[1],0]+ipos[nindex[1],0])
              elif N==2 and len(set(magnifications))==1:
                duplicated = True
                index = self.duplicate_locations(magnifications)
                ipos[index,0] = 0.0
                ipos[index,1] = 0.0
              if duplicated: 
                if vb: print 'reset image positions:',ipos
                magnifications = self.magnifications(ipos)
                parities = np.sign(magnifications)
                if vb: print 'reset image magnifications:',magnifications
          
          # Onwards!      
          if fail and keepgoing and humph < 42:
            humph += 1
          else:
            done = True # Hooray! Straight to successful return from here.

        return ipos,fail

    
# ----------------------------------------------------------------------------

    # NB: MUST BE SYNCHRONIZED WITH DEFLECTIONS() AND INVERSE_MAGNIFICATION_TENSORS()
    def potentials(self, imagepositions):
        ipos = np.atleast_2d(imagepositions)
        dpos = ipos - np.outer(np.ones(len(ipos)), self.position)
        r = np.sqrt(np.sum(dpos * dpos, axis=1))
        phi = np.arctan2(dpos[:,1], dpos[:,0])
        psis = self.einsteinradius * r + 0.5 * r * r * self.gamma * np.cos(2. * (phi - self.phi))
        return psis

# ----------------------------------------------------------------------------

    # input: image-plane positions shape (N, 2)
    # output: potential values at those positions shape (N, )
    def time_delays(self, sourceposition, imagepositions):
        N = len(imagepositions)
        dpos = np.outer(np.ones(N), sourceposition) - imagepositions
        return 0.5 * np.sum(dpos**2, axis=1) - self.potentials(imagepositions)

# ----------------------------------------------------------------------------

    # NB: MUST BE SYNCHRONIZED WITH POTENTIALS() AND INVERSE_MAGNIFICATION_TENSORS()
    def deflections(self, imagepositions):
        ipos = np.atleast_2d(imagepositions)
        dpos = ipos - np.outer(np.ones(len(ipos)), self.position)
        r = np.outer(np.sqrt(np.sum(dpos * dpos, axis=1)), [1,1])
        alphas = self.einsteinradius * dpos / r
        alphas[:,0] += self.gammacos2phi * dpos[:,0]
        alphas[:,0] += self.gammasin2phi * dpos[:,1]
        alphas[:,1] -= self.gammacos2phi * dpos[:,1]
        alphas[:,1] += self.gammasin2phi * dpos[:,0]
        return alphas

# ----------------------------------------------------------------------------

    # input: image positions shape (N, 2)
    # output: source position
    # note outer, sqrt, sum craziness
    def source_positions(self, imagepositions):
        return imagepositions - self.deflections(imagepositions)

 # ----------------------------------------------------------------------------

    # output shape (N, 2, 2)
    # NB: MUST BE SYNCHRONIZED WITH POTENTIALS() AND DEFLECTIONS()
    def inverse_magnification_tensors(self, imagepositions):
        tiny = 1e-12 # to deal with lens and image both exactly at origin
        ipos = np.atleast_2d(imagepositions)
        mag = np.zeros((len(ipos), 2, 2))
        mag[:,0,0] = 1.
        mag[:,1,1] = 1.
        # print "inverse_magnification_tensors: ipos = ",ipos
        dpos = ipos - self.position
        rcubed = np.sum(dpos * dpos, axis=1)**1.5 + tiny # tiny little hack
        if np.min(rcubed) <= 0.0:
            print "image positions = ",ipos
            print "lens position = ",self.position
            print "differences = ",dpos
            print "rcubed = ",rcubed
            print self
        else:
            mag[:,0,0] -= self.einsteinradius * dpos[:,1] * dpos[:,1] / rcubed
            mag[:,0,1] += self.einsteinradius * dpos[:,1] * dpos[:,0] / rcubed
            mag[:,1,0] += self.einsteinradius * dpos[:,0] * dpos[:,1] / rcubed
            mag[:,1,1] -= self.einsteinradius * dpos[:,0] * dpos[:,0] / rcubed
            mag[:,0,0] -= self.gammacos2phi
            mag[:,0,1] -= self.gammasin2phi
            mag[:,1,0] -= self.gammasin2phi
            mag[:,1,1] += self.gammacos2phi
        assert(np.min(rcubed) > 0.0)
        return mag

# ----------------------------------------------------------------------------

    def magnification_tensors(self, imagepositions):
        return np.array([np.linalg.inv(t) for t in self.inverse_magnification_tensors(imagepositions)])

    # crazy if you run this and magnificationtensors in the same code; think caching
    def magnifications(self, imagepositions):
        return 1. / np.array(map(np.linalg.det, self.inverse_magnification_tensors(imagepositions)))

    # crazy if you run this and magnificationtensors in the same code; think caching
    def parities(self, imagepositions):
        return np.sign(self.magnifications(imagepositions))

# ----------------------------------------------------------------------------

    # APPROXIMATION: Not yet marginalizing over true source position, true source flux
    # look for "WRONG" in code
    # Note magnification sign insanity
    def ln_prior(self, imagepositions, imagefluxes, positionvariance, fluxvariance, paritycheck=True):
        def ln_Gaussian_1d_zeromean(x, var):
            return -0.5 * np.log(2. * np.pi * var) - 0.5 * x**2 / var
        assert(len(imagepositions) == 4)
        sourcepositions = self.source_positions(imagepositions)
        meansourceposition = np.mean(sourcepositions, axis=0) # WRONG
        magtensors = self.magnification_tensors(imagepositions)
        mags = np.array(map(np.linalg.det, magtensors))
        if paritycheck:
            if mags[0] <= 0.:
                return -np.Inf
            if mags[1] >= 0.:
                return -np.Inf
            if mags[2] <= 0.:
                return -np.Inf
            if mags[3] >= 0.:
                return -np.Inf
        dimagepositions = np.array([np.dot(tens, (spos - meansourceposition)) for tens, spos in zip(magtensors, sourcepositions)])
        return np.sum(ln_Gaussian_1d_zeromean(dimagepositions, positionvariance))

# ----------------------------------------------------------------------------

    def sample_prior(self, imagepositions, imagefluxes, positionvariance, fluxvariance, nlink):
        def lnp(pars):
            return self.ln_prior(np.reshape(pars[:8], (4, 2)), None, positionvariance, None)
        pars = np.ravel(imagepositions)
        nwalkers = 100
        ndim     = len(pars)
        initial_position = [pars + 0.01 * np.random.normal(size=ndim) for i in xrange(nwalkers)]
        sampler = dfm.EnsembleSampler(nwalkers, ndim, lnp)
        pos, prob, state = sampler.run_mcmc(initial_position, None, nlink)
        print 'Mean acceptance fraction: ',np.mean(sampler.acceptance_fraction())
        return (sampler.get_chain(), sampler.get_lnprobability())

# ----------------------------------------------------------------------------

    # if timedelaymap then plot time-delay surface
    # if magnificationmap then plot magnification surface
    # several MAGIC NUMBERS
    def plot(self, sourcepositions=None, imagepositions=None, magnificationmap=False, timedelaymap=False):
        causticlw = 0.5
        tc = self.tangential_caustic().T
        plt.plot(tc[0], tc[1], 'k', lw=causticlw)
        rc = self.radial_caustic().T
        plt.plot(rc[0], rc[1], 'k', lw=causticlw)
        cc = self.critical_curve().T
        plt.plot(cc[0], cc[1], 'k', lw=3*causticlw)
        if sourcepositions is not None:
            spos = np.atleast_2d(sourcepositions)
            # print 'plot: plotting spos:', spos
            plt.scatter(spos[:,0], spos[:,1], c='k', marker='x', lw=causticlw)
        if imagepositions is not None:
            ipos = np.atleast_2d(imagepositions)
            mags = np.array(map(np.linalg.det, self.magnification_tensors(ipos)))
            s = 20. * np.sqrt(np.abs(mags))
            I = mags < 0
            if np.sum(I) > 0:
                plt.scatter(ipos[I,0], ipos[I,1], s=s[I], c='k', marker='o', facecolor='none')
            I = mags > 0
            if np.sum(I) > 0:
                plt.scatter(ipos[I,0], ipos[I,1], s=s[I], c='k', marker='s', facecolor='none')
        plt.xlabel('x (arcsec)')
        plt.ylabel('y (arcsec)')
        plt.axes().set_aspect('equal')
        plt.title('%s' % self)
        if magnificationmap:
            tcontours = np.arange(-100., 100.5, 1.0)
            xa, xb = plt.xlim()
            ya, yb = plt.ylim()
            xg, yg = np.meshgrid(np.arange(xa, xb, 0.01 * self.einsteinradius),
                                np.arange(ya, yb, 0.01 * self.einsteinradius))
            ig = np.array(zip(np.ravel(xg), np.ravel(yg)))
            tg = np.reshape(self.magnifications(ig), xg.shape)
            plt.contour(xg, yg, tg, tcontours, alpha=0.5, linewidths=causticlw)
        if timedelaymap:
            if sourcepositions is None:
                spos = np.atleast_2d(np.zeros(2))
            if imagepositions is None:
                ipos = self.guess_image_positions(spos)
            dts = self.time_delays(spos[0], ipos)
            ta = np.min(dts)
            tb = np.max(dts) + self.einsteinradius**2
            tcontours = np.arange(ta, tb, 0.01 * self.einsteinradius**2)
            xa, xb = plt.xlim()
            ya, yb = plt.ylim()
            xg, yg = np.meshgrid(np.arange(xa, xb, 0.01 * self.einsteinradius),
                                np.arange(ya, yb, 0.01 * self.einsteinradius))
            ig = np.array(zip(np.ravel(xg), np.ravel(yg)))
            tg = np.reshape(self.time_delays(spos[0], ig), xg.shape)
            plt.gray()
            plt.contour(xg, yg, tg, tcontours, alpha=0.5, linewidths=causticlw)
        return None

# ============================================================================

# TESTS:

# Draw a few hundred source positions close to a cusp, and see how the code 
# fails to find all images.

def merger_test(config):
    print "Requested image configuration is",config

    lenspos = [0.5, 0.75]
    b = 1.3 # arcsec
    if config == 'naked_cusp':
        gamma = 0.5
    else:
        gamma = 0.2
    phi = 0.2 # rad
    sis = GravitationalLens(lenspos, b, gamma, phi)

    plt.clf()
    sis.plot()
    foofile = '%s_diagram.png' % config
    plt.savefig(foofile)
    print "Lens outline plotted in",foofile

    nsample = 1000
    print "Drawing",nsample,"sample image positions..."
    ipos = np.zeros((nsample, 2))

    if config == 'minor_cusp':
        ipos[:,0] = np.random.uniform(1.5, 3.0, size=nsample)
        ipos[:,1] = np.random.uniform(0.5, 1.5, size=nsample)
    elif config == 'major_cusp':
        ipos[:,0] = np.random.uniform(0.5, 1.0, size=nsample)
        ipos[:,1] = np.random.uniform(-0.6, -1.0, size=nsample)
    elif config == 'naked_cusp':
        ipos[:,0] = np.random.uniform(0.5, 1.0, size=nsample)
        ipos[:,1] = np.random.uniform(-0.7, -2.0, size=nsample)
    else:
        ipos[:,0] = np.random.uniform(-2.0, 3.0, size=nsample)
        ipos[:,1] = np.random.uniform(-2.0, 3.0, size=nsample)
  
    # Mergers test - find an image position that does not solve properly, by
    # testing target and solved image numbers:
    spoz = sis.source_positions(ipos)
    nimz = np.ones(nsample)
    for i in range(nsample):
        nimz[i] = sis.number_of_images(spoz[i])
    if config == 'naked_cusp':
        target = 3
    else:
        target = 4
    index = np.where(nimz == target)
    print len(index[0]), 'randomly generated image positions hit target source area'
    ipoz = ipos[index]
    spoz = spoz[index]
    fail = False
    total = len(spoz)
    
    print 'Looping over these systems, looking for a failure:'
    keepgoing = False
    count = 0
    for spos in spoz:
        ipos,fail = sis.image_positions(spos,keepgoing)
        count += 1
        if fail:
          print "  Image solve failure (after ",int(count/(0.01*total)),"% of samples)"
          break
    
    if fail:
        print "  Source position is",spos
        nim = sis.number_of_images(spos)
        print "  No. of expected images:",nim
        ipos,dummy = sis.image_positions(spos,keepgoing)
        print "  Solved (but wrong) image positions:",ipos
        print "  (Wrong) image magnifications:",sis.magnifications(ipos)
        spos2 = sis.source_positions(ipos)
        print "  Corresponding source positions:",spos2
        plt.clf()
        sis.plot(sourcepositions=np.append(np.atleast_2d(spos), spos2, axis=0), imagepositions=ipos, timedelaymap=True)
        barfile = '%s_failure_%s.png' % (config,str(count))
        plt.savefig(barfile)
        print "  Problem images and source(s) plotted in",barfile
    
        keepgoing = True
        ipos,dummy = sis.image_positions(spos,keepgoing)
        print "  Solved (and correct) image positions:",ipos
        print "  (Correct) image magnifications:",sis.magnifications(ipos)
        spos2 = sis.source_positions(ipos)
        print "  Corresponding source positions:",spos2
        plt.clf()
        sis.plot(sourcepositions=np.append(np.atleast_2d(spos), spos2, axis=0), imagepositions=ipos, timedelaymap=True)
        barfile = '%s_success_%s.png' % (config,str(count))
        plt.savefig(barfile)
        print "  Problem images and source(s) plotted in",barfile
    
    else:
        print "All tests passed OK"

# ----------------------------------------------------------------------------

# Draw a few thousand source and lens configurations, and see how long the 
# code takes to find all images in each case.

import time

def speed_test():
    
    nsample = 1000
    print "Solving",nsample,"lens systems..."
    
    lenspos = [0.0, 0.0]
    b = 1.0 # arcsec

    gamma = np.random.uniform(0.0, 0.5, size=nsample)
    phi = np.random.uniform(0.0, np.pi, size=nsample)
    
    spos = np.zeros([nsample,2])
    spos[:,0] = np.random.uniform(-b, b, size=nsample)
    spos[:,1] = np.random.uniform(-b, b, size=nsample)
    
    nim = np.ones(nsample)
    dt = np.zeros(nsample)
    
    keepgoing = True
    nfail = 0
    for i in range(nsample):
        sis = GravitationalLens(lenspos, b, gamma[i], phi[i])
        nim[i] = sis.number_of_images(spos[i])
        if vb: print i, spos[i], gamma[i], phi[i], nim[i]
        t0 = time.time()
        ipos,fail = sis.image_positions(spos[i],keepgoing)
        dt[i] = time.time() - t0
        if vb: print ipos, dt[i]
        if fail:
            plt.clf()
            sis.plot(sourcepositions=spos[i],imagepositions=ipos,timedelaymap=True)
            foofile = 'tests/solving/failure%s_diagram.png' % str(i)
            plt.savefig(foofile)
            print "Solve failure: lens outline plotted in",foofile
            print "  image positions:",ipos
            print "  magnifications:",sis.magnifications(ipos)
            nfail += 1
        
    print "...done in %d seconds!" % sum(dt)
    plt.clf()
    plt.hist(dt, bins=10**np.linspace(-3, 1, num=20))
    plt.xscale('log')
    plt.xlabel('solve time (seconds)')
    plt.ylabel('frequency')
    histfile = 'tests/solving/speed_test_histogram.png'
    plt.savefig(histfile)
    print "Histogram plotted in",histfile
    
    print "Solve failure rate =",nfail*100.0/nsample,"%"


# ============================================================================

if __name__ == '__main__':

# Can we make the solver fail at merging image positions?
#     for config in ['major_cusp', 'minor_cusp', 'naked_cusp']:
#         merger_test(config)

# How often does solver fail, and how long does each solve take?
    speed_test()
    

# if False:
#     caustic = sis.tangential_caustic()
#     critcurve = sis.critical_curve()
#     posvar = 1.e-3
# 
#     imagefluxes = np.array([10., 10., 10., 10.])
#     magtensors = sis.magnificationtensors(imagepositions)
#     magscalars = np.array([np.linalg.det(ten) for ten in magtensors])
#     chain, lnp = sis.sample_prior(imagepositions, None, posvar, None, 10)
#     
#     plt.clf()
#     for i in range(5):
#         plt.plot(lnp[i,:])
#     plt.xlabel('link number')
#     plt.ylabel('ln prior probability')
#     plt.savefig('lnp.png')
#     for i in range(5):
#         if i == 0:
#             ipos = imagepositions
#             print sis.ln_prior(imagepositions, None, posvar, None)
#         else:
#             ipos = np.reshape(chain[i,:,-1], (4,2))
#         spos = sis.sourcepositions(ipos)
#         plt.clf()
#         plt.subplot(111, aspect='equal')
#         mags = np.abs(map(np.linalg.det, sis.magnificationtensors(ipos)))
#         plt.plot(caustic[:,0], caustic[:,1], 'g')
#         plt.plot(critcurve[:,0], critcurve[:,1], 'b')
#         plt.scatter(ipos[:,0], ipos[:,1], s=mags, c='k', marker='o')
#         plt.scatter(spos[:,0], spos[:,1], c='r', marker='o')
#         xt = 3 * np.random.uniform(size=(3000,2))
#         color = ['m', 'r', 'k', 'g', 'c', 'b']
#         colors = [color[sis.number_of_images(x)] for x in xt]
#         plt.scatter(xt[:,0], xt[:,1], c=colors, marker='o', alpha=0.5)
#         plt.xlim(-2, 3.)
#         plt.ylim(-2, 3.)
#         plt.xlabel('x (arcsec)')
#         plt.ylabel('y (arcsec)')
#         plt.savefig('sky-%02d.png' % i)
# 
# ============================================================================
