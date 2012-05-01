'''
This file is part of the LensFinder project.
Copyright 2011 2012 David W. Hogg (NYU) and Phil Marshall (Oxford).

Description
-----------

Nebula is a flexible model for a generic extended object, chosen to be the sum of np 
point source components and ne extended source components (composite galaxies)


To-do
-----
- inherit tractor star/galaxy properties
- enable np, ne sources... where np is 4 and ne = 1.
'''

import numpy as np
# import markovpy as dfm
if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import rc
    rc('font',**{'family':'serif','serif':'Computer Modern Roman','size':18})
    rc('text', usetex=True)
    import pylab as plt
import matplotlib.nxutils as nx

# ============================================================================

# flexible extended, multicomponent model
# initialize with
# - position        : galaxy position (same units as extent)
# - np              : number of point sources (which can be variable)
# - ne              : number of extended sources (which cannot)
# - extent          : radius within which all sources are centred (some angular units)

class Nebula:

    def __init__(self, position, extent, np, ne):
        self.name = 'Nebula'
        self.position = np.array(position)
        self.extent = extent
        self.np = np
        self.ne = ne
        return None

    def __str__(self):
        return '%s, position (%.2f,%.2f), extent %.2f, Np,Ne = %d, %d)' % (self.name, self.position, self.extent, self.np, self.ne)

# ----------------------------------------------------------------------------

    # Accept source, each one a Gaussian
    def image_positions(self, sourceposition):
        
        return None

# ----------------------------------------------------------------------------
# 
#     # NB: MUST BE SYNCHRONIZED WITH DEFLECTIONS() AND INVERSE_MAGNIFICATION_TENSORS()
#     def potentials(self, imagepositions):
#         ipos = np.atleast_2d(imagepositions)
#         dpos = ipos - np.outer(np.ones(len(ipos)), self.position)
#         r = np.sqrt(np.sum(dpos * dpos, axis=1))
#         phi = np.arctan2(dpos[:,1], dpos[:,0])
#         psis = self.einsteinradius * r + 0.5 * r * r * self.gamma * np.cos(2. * (phi - self.phi))
#         return psis
# 
#     # input: image-plane positions shape (N, 2)
#     # output: potential values at those positions shape (N, )
#     def time_delays(self, sourceposition, imagepositions):
#         N = len(imagepositions)
#         dpos = np.outer(np.ones(N), sourceposition) - imagepositions
#         return 0.5 * np.sum(dpos**2, axis=1) - self.potentials(imagepositions)
# 
# # ----------------------------------------------------------------------------
# 
#     # NB: MUST BE SYNCHRONIZED WITH POTENTIALS() AND INVERSE_MAGNIFICATION_TENSORS()
#     def deflections(self, imagepositions):
#         ipos = np.atleast_2d(imagepositions)
#         dpos = ipos - np.outer(np.ones(len(ipos)), self.position)
#         r = np.outer(np.sqrt(np.sum(dpos * dpos, axis=1)), [1,1])
#         alphas = self.einsteinradius * dpos / r
#         alphas[:,0] += self.gammacos2phi * dpos[:,0]
#         alphas[:,0] += self.gammasin2phi * dpos[:,1]
#         alphas[:,1] -= self.gammacos2phi * dpos[:,1]
#         alphas[:,1] += self.gammasin2phi * dpos[:,0]
#         return alphas
# 
#     # input: image positions shape (N, 2)
#     # output: source position
#     # note outer, sqrt, sum craziness
#     def source_positions(self, imagepositions):
#         return imagepositions - self.deflections(imagepositions)
# 
#  # ----------------------------------------------------------------------------
# 
#     # output shape (N, 2, 2)
#     # NB: MUST BE SYNCHRONIZED WITH POTENTIALS() AND DEFLECTIONS()
#     def inverse_magnification_tensors(self, imagepositions):
#         ipos = np.atleast_2d(imagepositions)
#         mag = np.zeros((len(ipos), 2, 2))
#         mag[:,0,0] = 1.
#         mag[:,1,1] = 1.
#         # print "inverse_magnification_tensors: ipos = ",ipos
#         dpos = ipos - self.position
#         rcubed = np.sum(dpos * dpos, axis=1)**1.5
#         if np.min(rcubed) <= 0.:
#             print ipos
#             print self.position
#             print mag
#             print self
#         assert(np.min(rcubed) > 0.)
#         mag[:,0,0] -= self.einsteinradius * dpos[:,1] * dpos[:,1] / rcubed
#         mag[:,0,1] += self.einsteinradius * dpos[:,1] * dpos[:,0] / rcubed
#         mag[:,1,0] += self.einsteinradius * dpos[:,0] * dpos[:,1] / rcubed
#         mag[:,1,1] -= self.einsteinradius * dpos[:,0] * dpos[:,0] / rcubed
#         mag[:,0,0] -= self.gammacos2phi
#         mag[:,0,1] -= self.gammasin2phi
#         mag[:,1,0] -= self.gammasin2phi
#         mag[:,1,1] += self.gammacos2phi
#         return mag
# 
#     def magnification_tensors(self, imagepositions):
#         return np.array([np.linalg.inv(t) for t in self.inverse_magnification_tensors(imagepositions)])
# 
#     # crazy if you run this and magnificationtensors in the same code; think caching
#     def magnifications(self, imagepositions):
#         return 1. / np.array(map(np.linalg.det, self.inverse_magnification_tensors(imagepositions)))
# 
#     # crazy if you run this and magnificationtensors in the same code; think caching
#     def parities(self, imagepositions):
#         return np.sign(self.magnifications(imagepositions))
# 
# # ----------------------------------------------------------------------------
# 
#     # APPROXIMATION: Not yet marginalizing over true source position, true source flux
#     # look for "WRONG" in code
#     # Note magnification sign insanity
#     def ln_prior(self, imagepositions, imagefluxes, positionvariance, fluxvariance, paritycheck=True):
#         def ln_Gaussian_1d_zeromean(x, var):
#             return -0.5 * np.log(2. * np.pi * var) - 0.5 * x**2 / var
#         assert(len(imagepositions) == 4)
#         sourcepositions = self.source_positions(imagepositions)
#         meansourceposition = np.mean(sourcepositions, axis=0) # WRONG
#         magtensors = self.magnification_tensors(imagepositions)
#         mags = np.array(map(np.linalg.det, magtensors))
#         if paritycheck:
#             if mags[0] <= 0.:
#                 return -np.Inf
#             if mags[1] >= 0.:
#                 return -np.Inf
#             if mags[2] <= 0.:
#                 return -np.Inf
#             if mags[3] >= 0.:
#                 return -np.Inf
#         dimagepositions = np.array([np.dot(tens, (spos - meansourceposition)) for tens, spos in zip(magtensors, sourcepositions)])
#         return np.sum(ln_Gaussian_1d_zeromean(dimagepositions, positionvariance))
# 
# # ----------------------------------------------------------------------------
# 
#     def sample_prior(self, imagepositions, imagefluxes, positionvariance, fluxvariance, nlink):
#         def lnp(pars):
#             return self.ln_prior(np.reshape(pars[:8], (4, 2)), None, positionvariance, None)
#         pars = np.ravel(imagepositions)
#         nwalkers = 100
#         ndim     = len(pars)
#         initial_position = [pars + 0.01 * np.random.normal(size=ndim) for i in xrange(nwalkers)]
#         sampler = dfm.EnsembleSampler(nwalkers, ndim, lnp)
#         pos, prob, state = sampler.run_mcmc(initial_position, None, nlink)
#         print 'Mean acceptance fraction: ',np.mean(sampler.acceptance_fraction())
#         return (sampler.get_chain(), sampler.get_lnprobability())
# 
# # ----------------------------------------------------------------------------
# 
#     # if timedelaymap then plot time-delay surface
#     # if magnificationmap then plot magnification surface
#     # several MAGIC NUMBERS
#     def plot(self, sourcepositions=None, imagepositions=None, magnificationmap=False, timedelaymap=False):
#         causticlw = 0.5
#         tc = self.tangential_caustic().T
#         plt.plot(tc[0], tc[1], 'k', lw=causticlw)
#         rc = self.radial_caustic().T
#         plt.plot(rc[0], rc[1], 'k', lw=causticlw)
#         cc = self.critical_curve().T
#         plt.plot(cc[0], cc[1], 'k', lw=3*causticlw)
#         if sourcepositions is not None:
#             spos = np.atleast_2d(sourcepositions)
#             # print 'plot: plotting spos:', spos
#             plt.scatter(spos[:,0], spos[:,1], c='k', marker='x', lw=causticlw)
#         if imagepositions is not None:
#             ipos = np.atleast_2d(imagepositions)
#             mags = np.array(map(np.linalg.det, self.magnification_tensors(ipos)))
#             s = 20. * np.sqrt(np.abs(mags))
#             I = mags < 0
#             if np.sum(I) > 0:
#                 plt.scatter(ipos[I,0], ipos[I,1], s=s[I], c='k', marker='o', facecolor='none')
#             I = mags > 0
#             if np.sum(I) > 0:
#                 plt.scatter(ipos[I,0], ipos[I,1], s=s[I], c='k', marker='s', facecolor='none')
#         plt.xlabel('x (arcsec)')
#         plt.ylabel('y (arcsec)')
#         plt.axes().set_aspect('equal')
#         plt.title('%s' % self)
#         if magnificationmap:
#             tcontours = np.arange(-100., 100.5, 1.0)
#             xa, xb = plt.xlim()
#             ya, yb = plt.ylim()
#             xg, yg = np.meshgrid(np.arange(xa, xb, 0.01 * self.einsteinradius),
#                                 np.arange(ya, yb, 0.01 * self.einsteinradius))
#             ig = np.array(zip(np.ravel(xg), np.ravel(yg)))
#             tg = np.reshape(self.magnifications(ig), xg.shape)
#             plt.contour(xg, yg, tg, tcontours, alpha=0.5, linewidths=causticlw)
#         if timedelaymap:
#             dts = self.time_delays(spos[0], ipos)
#             ta = np.min(dts)
#             tb = np.max(dts) + self.einsteinradius**2
#             tcontours = np.arange(ta, tb, 0.01 * self.einsteinradius**2)
#             xa, xb = plt.xlim()
#             ya, yb = plt.ylim()
#             xg, yg = np.meshgrid(np.arange(xa, xb, 0.01 * self.einsteinradius),
#                                 np.arange(ya, yb, 0.01 * self.einsteinradius))
#             ig = np.array(zip(np.ravel(xg), np.ravel(yg)))
#             tg = np.reshape(self.time_delays(spos[0], ig), xg.shape)
#             plt.gray()
#             plt.contour(xg, yg, tg, tcontours, alpha=0.5, linewidths=causticlw)
#         return None
# 
# ============================================================================

# 

def test():
    
    print "Generating a random nebula:"


       light = stgal.CompositeGalaxy(pos, brexp, shexp, brdev, shdev)
       mass = ....
       lens = LensGalaxy(light, mass)
       source = PointSource(pos, bright)

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
    foofile = 'foo_%s.png' % config
    plt.savefig(foofile)
    print "Lens outline plotted in",foofile

    nsample = 100
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
    print len(index[0]), 'randomly generated image positions hit target'
    ipoz = ipos[index]
    spoz = spoz[index]
    fail = False
    total = len(spoz)
    count = 0
    for spos in spoz:
        ipos,fail = sis.image_positions(spos)
        count += 1
        if fail:
            print "Image solve failure (after ",int(count/(0.01*total)),"% of samples)"
            break
    if fail:
        print "Source position is",spos
        nim = sis.number_of_images(spos)
        print "No. of expected images:",nim
        ipos,dummy = sis.image_positions(spos)
        print "Solved image positions:",ipos
        print "Image magnifications:",sis.magnifications(ipos)
        spos2 = sis.source_positions(ipos)
        print "Corresponding source positions:",spos2
        plt.clf()
        sis.plot(sourcepositions=np.append(np.atleast_2d(spos), spos2, axis=0), imagepositions=ipos, timedelaymap=True)
        barfile = 'bar_%s.png' % config
        plt.savefig(barfile)
        print "Images and source(s) plotted in",barfile
    else:
        print "All tests passed OK"

# ============================================================================

if __name__ == '__main__':
    
    test()

# ============================================================================
