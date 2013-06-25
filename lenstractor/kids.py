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

def KIDSWCS(hdr):
      '''
      Return a WCS object initialised from a KIDS file header.
      WARNING: PC matrix not being used, determinant may be wrong sign...
      Need to check with literature image of H1413

      TODO: CD1_1 has a negative value in KIDS where CDELT1 in PS1 is
            positive. Need to check whether this is correct.
      '''
      t = util.Tan()
      t.set_crpix(hdr['CRPIX1'], hdr['CRPIX2'])
      t.set_crval(hdr['CRVAL1'], hdr['CRVAL2'])
      t.set_cd(hdr['CD1_1'], 0., 0., hdr['CD2_2'])
      t.set_imagesize(hdr['NAXIS1'], hdr['NAXIS2'])
      
      return tractor.FitsWcs(t)


# ============================================================================

def KIDS_IQ(hdr):
   """
   Return the image quality (FWHM) in pixels given a FITS file header.
   """
   # TODO: Determine the correct units of PSF_RAD.
   FWHM = hdr['PSF_RAD'] # arcsec
   # Need it in pixels:
   # TODO: get the pixel scale from the headers
   FWHM /= 0.2
   
   return FWHM

# ============================================================================

def KIDS_photocal(hdr):
   """
   Return a PS1MagsPhotoCal object given a FITS file header.
   """
   # TODO: The KiDS zeropoint is always 0 in CoaddedRegriddedFrames.
   #       This might not give the desired result.
   band = hdr['FILT_ID'][5] # e.g. 'OCAM_g_SDSS'[5] -> 'g'
   zpt = hdr['ZEROPNT']
   photocal = lenstractor.PS1MagsPhotoCal(zpt,band)

   return band,photocal


# ============================================================================

def KIDS_epoch(hdr):
   """
   Return the MJD given a FITS file header.
   """
   # TODO: OBS_STRT has a value like '2013-04-07T06:21:19' in KiDS, while
   #       the 'MJD-OBS' in PS1 is a number converted to a string.
   return hdr['OBS_STRT']

