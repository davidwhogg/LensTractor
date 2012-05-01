'''
This file is part of the LensFinder project.
Copyright 2012 David W. Hogg (NYU) and Phil Marshall (Oxford).

Description
-----------

Testing tractor functionality - read in an image and its weight map,
guess the PSF, and make a point source image on the same grid and take
the difference.

'''

import numpy as np, os

import tractor

import lensfinder

# ============================================================================

def main():
   
   # Read in a deck of postcard images:
   
   folder = os.environ['LENSFINDER_DIR']+'/examples'
   postcards = lensfinder.deck(folder)
   
   
#    W,H = 300,200
# 
#    psf = NCircularGaussianPSF([2.], [1.])
#    cx,cy,flux = 100., 80., 1000.
# 
#    image = np.zeros((H,W))
#    err = np.zeros_like(image) + 1.
#    invvar = 1./(err**2)
#    src = PointSource(PixPos(cx, cy), Flux(flux))
#    photocal = NullPhotoCal()
#    wcs = NullWCS()
# 
#    data = Image(data=image, invvar=invvar,
#                       psf=psf, wcs=wcs, photocal=photocal)
# 
#    # add perfect point source to image.
#    patch = src.getModelPatch(data)
#    patch.addTo(image)
# 
#    assert(abs(image.sum() - flux) < 1.)
# 
#    # Create new Image with the synthetic image.
#    data = Image(data=image, invvar=invvar,
#                       psf=psf, wcs=wcs, photocal=photocal)

# ============================================================================

if __name__ == '__main__':
   main()
