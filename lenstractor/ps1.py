# ============================================================================

'''
This file is part of the lenstractor project.
Copyright 2012 David W. Hogg (NYU) and Phil Marshall (Oxford).

Description
-----------

Data management classes and functions for the PS1 survey.

* Read in a deck of postcard images in FITS files, match them up and return an 
  array of tractor image data structures. 
* Manipulate PS1 FITS image headers, WCS and photo cal items.
'''

import numpy as np
import os,glob,string

from astrometry.util import util
import tractor
import lenstractor

vb = True

# ============================================================================

def PS1WCS(hdr):
      '''
      Return a WCS object initialised from a PS1 file header.
      WARNING: PC matrix not being used, determinant may be wrong sign...
      Need to check with literature image of H1413
      '''
      t = util.Tan()
      t.set_crpix(hdr['CRPIX1'], hdr['CRPIX2'])
      t.set_crval(hdr['CRVAL1'], hdr['CRVAL2'])
      t.set_cd(hdr['CDELT1'], 0., 0., hdr['CDELT2'])
      t.set_imagesize(hdr['NAXIS1'], hdr['NAXIS2'])
      
      return tractor.FitsWcs(t)
         
# ============================================================================

class PS1MagsPhotoCal(tractor.BaseParams):
      '''
      A photocal for a Mags brightness object.
      '''
      def __init__(self, zpt, bandname):
            self.bandname = bandname
            self.zpt = zpt

      def __str__(self):
            return (self.getName()+': '+self.bandname+' band, zpt='+str(self.getParams()))
      
      @staticmethod
      def getName():
            return 'PS1MagsPhotoCal'
      def getParams(self):
            return [self.zpt]
      def getStepSizes(self):
            return [0.01]
      def setParam(self, i, p):
            assert(i == 0)
            self.zpt = p
      def getParamNames(self):
            return ['zpt']
      def hashkey(self):
            return (self.getName(), self.bandname, self.zpt)

      def brightnessToCounts(self, brightness):
            mag = brightness.getMag(self.bandname)
            if not np.isfinite(mag):
                  return 0.
            # MAGIC
            if mag > 50.:
                  return 0.
            # Assume zeropoint is the whole story:
            return 10.**(-0.4 * (mag - self.zpt))

# ============================================================================
# Return a PS1MagsPhotoCal object given a FITS file header.

def PS1_photocal(hdr):

   band = hdr['HIERARCH FPA.FILTER'][0]
   zpt = hdr['HIERARCH FPA.ZP']
   photocal = lenstractor.PS1MagsPhotoCal(zpt,band)

   return band,photocal

# ============================================================================
# Return the MJD given a FITS file header.

def PS1_epoch(hdr):

   return "%.5f" % hdr['MJD-OBS']

# ============================================================================
# Return the image quality (FWHM) in pixels given a FITS file header.

def PS1_IQ(hdr):

   FWHM = hdr['HIERARCH CHIP.SEEING']
   if FWHM == 'NaN': 
      FWHM = 1.0 # arcsec
      # Need it in pixels:
      FWHM = FWHM/(3600.0*hdr['CDELT1'])
      print "PS1_IQ: WARNING: FWHM = NaN in header, setting to 1.0 arcsec =",FWHM,"pixels"
   
   
   return FWHM

# ============================================================================

if __name__ == '__main__':

     if False:
        # Testing LensPlaneWCS:

        ra,dec = 310.0,0.0
        pos = tractor.RaDecPos(ra,dec)

        lenswcs = lenstractor.LensPlaneWCS(pos)

        print lenswcs

        dt = 1/360.0
        print lenswcs.positionToPixel(tractor.RaDecPos(ra, dec))
        print lenswcs.positionToPixel(tractor.RaDecPos(ra+dt, dec))
        print lenswcs.positionToPixel(tractor.RaDecPos(ra-dt, dec))
        print lenswcs.positionToPixel(tractor.RaDecPos(ra, dec+dt))
        print lenswcs.positionToPixel(tractor.RaDecPos(ra, dec-dt))

      
# ============================================================================
