'''
This file is part of the LensFinder project.
Copyright 2012 David W. Hogg (NYU) and Phil Marshall (Oxford).

Description
-----------

Running the tractor on PS1 *single object cutout* images.
Read in an image and its weight map, guess the PSF, put an object at the
centre image of the field and update the catalog.

Example use
-----------

 python ps1tractor.py -v \
    examples/H1413+117_10x10arcsec_55377.34051_z_sci.fits \
    examples/H1413+117_10x10arcsec_55377.34051_z_var.fits

'''

if __name__ == '__main__':
      import matplotlib
      matplotlib.use('Agg')

import os
import logging
import numpy as np
import pylab as plt
import pyfits

from astrometry.util.file import *
# from astrometry.util.util import Sip
from astrometry.util.pyfits_utils import *

import tractor
from tractor import sdss_galaxy as gal
import lensfinder

# ============================================================================

def ps1tractor():

   from optparse import OptionParser
   import sys

   # Set available options:
   parser = OptionParser(usage=('%prog <sci> <var>'))
   # Verbosity:
   parser.add_option('-v', '--verbose', dest='verbose', action='count', \
                     default=False, help='Make more verbose')
   
   # Read in options and arguments - note only sci and wht images are supplied:
   opt,args = parser.parse_args()
   
   if len(args) != 2:
      parser.print_help()
      sys.exit(-1)
   scifile, varfile = args
 
   # -------------------------------------------------------------------------
   # Logging to terminal:
   
   if opt.verbose:
      lvl = logging.DEBUG
   else:
      lvl = logging.INFO
   logging.basicConfig(level=lvl, format='%(message)s', stream=sys.stdout)

   # -------------------------------------------------------------------------
   # Get sci and wht images, and make mask:
   
   hdulist = pyfits.open(scifile)
   sci = hdulist[0].data
   hdr = hdulist[0].header
   hdulist.close()
   NX,NY = sci.shape

   hdulist = pyfits.open(varfile)
   var = hdulist[0].data
   hdulist.close()
   assert sci.shape == var.shape
   
   mask = numpy.ones([NX,NY],dtype=np.int16)
   mask[numpy.where(var == 0)] = 0
   
   # Assign var=nan all zero weight:
   var[var != var] = 0.0
   # Convert var to wht, and find median uncertainty as well:
   invvar = 1.0/var
   # Assign var=0, var<0 all zero weight:
   invvar[var == 0] = 0.0
   invvar[var < 0] = 0.0
   # Zero out sci image where wht is 0.0:
   sci[invvar == 0] = 0.0
   # Rough estimate of global pixel uncertainty:
   sig = np.sqrt(np.median(var))
   
   # Report on progress so far:
   if opt.verbose:
     print 'Sci header:', hdr
     print 'Read in sci image:', sci.shape, sci
     print 'Read in var image:', var.shape, var
     print 'Made mask image:', mask.shape, mask
     print 'Variance range:', var.min(), var.max()
     print 'Median pixel uncertainty:', sig
     print 'Image median:', np.median(sci.ravel())

   # -------------------------------------------------------------------------

#    for k,v in maskplanes.items():
#          plt.clf()
# 
#          I = ((mask & (1 << v)) != 0)
#          rgb = np.zeros((NX,NY,3))
#          clipimg = np.clip((img - (-3.*sig)) / (13.*sig), 0, 1)
#          cimg = clipimg.copy()
#          cimg[I] = 1
#          rgb[:,:,0] = cimg
#          cimg = clipimg.copy()
#          cimg[I] = 0
#          rgb[:,:,1] = cimg
#          rgb[:,:,2] = cimg
#          plt.imshow(rgb, interpolation='nearest', origin='lower')
#          plt.title(k)
#          plt.savefig('mask-%s.png' % k.lower())
# 
#    badmask = sum([(1 << maskplanes[k]) for k in ['BAD', 'SAT', 'INTRP', 'CR']])
#    # HACK -- left EDGE sucks
#    badmask += (1 << maskplanes['EDGE'])
#    #badmask = (1 << 0) | (1 << 1) | (1 << 2) | (1 << 3)
#    #badmask |= (1 << 4)
#    print 'Masking out: 0x%x' % badmask
#    invvar[(mask & badmask) != 0] = 0.
# 
#    assert(all(np.isfinite(img.ravel())))
#    assert(all(np.isfinite(invvar.ravel())))
# 
#    psf = pyfits.open(psffn)[0].data
#    print 'psf', psf.shape
#    psf /= psf.sum()
# 
#    from tractor.emfit import em_fit_2d
#    from tractor.fitpsf import em_init_params
# 
#    # Create Gaussian mixture model PSF approximation.
#    S = psf.shape[0]
#    # number of Gaussian components
#    K = 3
#    w,mu,sig = em_init_params(K, None, None, None)
#    II = psf.copy()
#    II /= II.sum()
#    # HIDEOUS HACK
#    II = np.maximum(II, 0)
#    print 'Multi-Gaussian PSF fit...'
#    xm,ym = -(S/2), -(S/2)
#    em_fit_2d(II, xm, ym, w, mu, sig)
#    print 'w,mu,sig', w,mu,sig
#    mypsf = tractor.GaussianMixturePSF(w, mu, sig)
# 
# 
#    P = mypsf.getPointSourcePatch(S/2, S/2)
#    mn,mx = psf.min(), psf.max()
#    ima = dict(interpolation='nearest', origin='lower',
#                   vmin=mn, vmax=mx)
#    plt.clf()
#    plt.subplot(1,2,1)
#    plt.imshow(psf, **ima)
#    plt.subplot(1,2,2)
#    pimg = np.zeros_like(psf)
#    P.addTo(pimg)
#    plt.imshow(pimg, **ima)
#    plt.savefig('psf.png')
# 
#    sig = np.sqrt(np.median(var))
# 
#    plt.clf()
#    plt.hist(img.ravel(), 100, range=(-3.*sig, 3.*sig))
#    plt.savefig('imghist.png')
# 
#    srcs = fits_table(srcfn)
#    print 'Initial:', len(srcs), 'sources'
#    # Trim sources with x=0 or y=0
#    srcs = srcs[(srcs.x != 0) * (srcs.y != 0)]
#    print 'Trim on x,y:', len(srcs), 'sources left'
#    # Zero out nans & infs
#    for c in ['theta', 'a', 'b']:
#          I = np.logical_not(np.isfinite(srcs.get(c)))
#          srcs.get(c)[I] = 0.
#    # Set sources with flux=NaN to something more sensible...
#    I = np.logical_not(np.isfinite(srcs.flux))
#    srcs.flux[I] = 1.
#    # Sort sources by flux.
#    srcs = srcs[np.argsort(-srcs.flux)]
# 
#    # Trim sources that are way outside the image.
#    margin = 8. * np.maximum(srcs.a, srcs.b)
#    H,W = img.shape
#    srcs = srcs[(srcs.x > -margin) * (srcs.y > -margin) *
#                      (srcs.x < (W+margin) * (srcs.y < (H+margin)))]
#    print 'Trim out-of-bounds:', len(srcs), 'sources left'
# 
# 
#    wcs = tractor.FitsWcs(Sip(imgfn, 1))
#    #wcs = tractor.NullWCS()
# 
#    timg = tractor.Image(data=img, invvar=invvar, psf=mypsf, wcs=wcs,
#                                   sky=tractor.ConstantSky(0.),
#                                   photocal=tractor.NullPhotoCal(),
#                                   name='image')
# 
#    inverr = timg.getInvError()
#    assert(all(np.isfinite(inverr.ravel())))
# 
#    tsrcs = []
#    for s in srcs:
#          #pos = tractor.PixPos(s.x, s.y)
#          pos = tractor.RaDecPos(s.ra, s.dec)
#          if s.a > 0 and s.b > 0:
#                eflux = tractor.Flux(s.flux / 2.)
#                dflux = tractor.Flux(s.flux / 2.)
#                re,ab,phi = s.a, s.b/s.a, 90.-s.theta
#                eshape = gal.GalaxyShape(re,ab,phi)
#                dshape = gal.GalaxyShape(re,ab,phi)
#                print 'Fluxes', eflux, dflux
#                tsrc = gal.CompositeGalaxy(pos, eflux, eshape, dflux, dshape)
#          else:
#                flux = tractor.Flux(s.flux)
#                print 'Flux', flux
#                tsrc = tractor.PointSource(pos, flux)
#          tsrcs.append(tsrc)
# 
#    chug = tractor.Tractor([timg])
#    for src in tsrcs:
#          if chug.getModelPatch(timg, src) is None:
#                print 'Dropping non-overlapping source:', src
#                continue
#          chug.addSource(src)
#    print 'Kept a total of', len(chug.catalog), 'sources'
# 
#    ima = dict(interpolation='nearest', origin='lower',
#                   vmin=-3.*sig, vmax=10.*sig)
#    chia = dict(interpolation='nearest', origin='lower',
#                      vmin=-5., vmax=5.)
# 
#    plt.clf()
#    plt.imshow(img, **ima)
#    plt.colorbar()
#    plt.savefig('img.png')
# 
#    plt.clf()
#    plt.imshow(invvar, interpolation='nearest', origin='lower')
#    plt.colorbar()
#    plt.savefig('invvar.png')
# 
#    mod = chug.getModelImages()[0]
#    plt.clf()
#    plt.imshow(mod, **ima)
#    plt.colorbar()
#    plt.savefig('mod-0.png')
# 
#    chi = chug.getChiImage(0)
#    plt.clf()
#    plt.imshow(chi, **chia)
#    plt.colorbar()
#    plt.savefig('chi-0.png')
# 
#    for step in range(5):
#          cat = chug.getCatalog()
#          for src in cat:
#                if chug.getModelPatch(timg, src) is None:
#                      print 'Dropping non-overlapping source:', src
#                      chug.removeSource(src)
#          print 'Kept a total of', len(chug.catalog), 'sources'
# 
#          #cat = chug.getCatalog()
#          #for i,src in enumerate([]):
#          #for i,src in enumerate(chug.getCatalog()):
#          #for i in range(len(cat)):
#          i = 0
#          while i < len(cat):
#                src = cat[i]
# 
#                #print 'Step', i
#                #for j,s in enumerate(cat):
#                #     x,y = timg.getWcs().positionToPixel(s, s.getPosition())
#                #     print '  ',
#                #     if j == i:
#                #           print '*',
#                #     print '(%6.1f, %6.1f)'%(x,y), s
# 
#                print 'Optimizing source', i, 'of', len(cat)
# 
#                x,y = timg.getWcs().positionToPixel(src, src.getPosition())
#                print '(%6.1f, %6.1f)'%(x,y), src
#                # pre = src.getModelPatch(timg)
# 
#                s1 = str(src)
#                print 'src1 ', s1
#                dlnp1,X,a = chug.optimizeCatalogFluxes(srcs=[src])
#                s2 = str(src)
#                dlnp2,X,a = chug.optimizeCatalogAtFixedComplexityStep(srcs=[src], sky=False)
#                s3 = str(src)
# 
#                #post = src.getModelPatch(timg)
# 
#                print 'src1 ', s1
#                print 'src2 ', s2
#                print 'src3 ', s3
#                print 'dlnp', dlnp1, dlnp2
# 
#                if chug.getModelPatch(timg, src) is None:
#                      print 'After optimizing, no overlap!'
#                      print 'Removing source', src
#                      chug.removeSource(src)
#                      i -= 1
#                i += 1
# 
#                # plt.clf()
#                # plt.subplot(2,2,1)
#                # img = timg.getImage()
#                # (x0,x1,y0,y1) = pre.getExtent()
#                # plt.imshow(img, **ima)
#                # ax = plt.axis()
#                # plt.plot([x0,x0,x1,x1,x0], [y0,y1,y1,y0,y0], 'k-', lw=2)
#                # plt.axis(ax)
#                # plt.subplot(2,2,3)
#                # plt.imshow(pre.getImage(), **ima)
#                # plt.subplot(2,2,4)
#                # plt.imshow(post.getImage(), **ima)
#                # plt.savefig('prepost-s%i-s%03i.png' % (step, i))
#                #
#                # mod = chug.getModelImages()[0]
#                # plt.clf()
#                # plt.imshow(mod, **ima)
#                # plt.colorbar()
#                # plt.savefig('mod-s%i-s%03i.png' % (step, i))
#                # chi = chug.getChiImage(0)
#                # plt.clf()
#                # plt.imshow(chi, **chia)
#                # plt.colorbar()
#                # plt.savefig('chi-s%i-s%03i.png' % (step, i))
# 
# 
#          #dlnp,x,a = chug.optimizeCatalogFluxes()
#          #print 'fluxes: dlnp', dlnp
#          #dlnp,x,a = chug.optimizeCatalogAtFixedComplexityStep()
#          #print 'opt: dlnp', dlnp
# 
#          mod = chug.getModelImages()[0]
#          plt.clf()
#          plt.imshow(mod, **ima)
#          plt.colorbar()
#          plt.savefig('mod-%i.png' % (step+1))
# 
#          chi = chug.getChiImage(0)
#          plt.clf()
#          plt.imshow(chi, **chia)
#          plt.colorbar()
#          plt.savefig('chi-%i.png' % (step+1))
# 
#    return
# 
#    for step in range(5):
#          chug.optimizeCatalogFluxes()
#          mod = chug.getModelImages()[0]
#          plt.clf()
#          plt.imshow(mod, **ima)
#          plt.colorbar()
#          plt.savefig('mod-s%i.png' % step)
# 
#          chi = chug.getChiImage(0)
#          plt.clf()
#          plt.imshow(chi, **chia)
#          plt.colorbar()
#          plt.savefig('chi-s%i.png' % step)
# 
# ============================================================================

if __name__ == '__main__':
   
   ps1tractor()

