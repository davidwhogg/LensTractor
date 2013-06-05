# ============================================================================

'''
This file is part of the lenstractor project.
Copyright 2012 David W. Hogg (NYU) and Phil Marshall (Oxford).
Copyright 2013 Hugo Buddelmeijer & Casper Blokzijl.

Description
-----------

Data management classes and functions for the KIDS survey.

* Read in a deck of postcard images in FITS files, match them up and return an 
  array of tractor image data structures. 
* Manipulate KIDS FITS image headers, WCS and photo cal items.
'''

import numpy as np
import os,glob,string

from astrometry.util import util
import tractor
import lenstractor

vb = True

# ============================================================================

# ============================================================================
# Return the image quality (FWHM) in pixels given a FITS file header.

def KIDS_IQ(hdr):
   # TODO: Determine the correct units of PSF_RAD.
   FWHM = hdr['PSF_RAD'] # arcsec
   # Need it in pixels:
   # TODO: get the pixel scale from the headers
   FWHM /= 0.2
   
   return FWHM
