# ============================================================================

'''
This file is part of the LensFinder project.
Copyright 2012 David W. Hogg (NYU) and Phil Marshall (Oxford).

Description
-----------

Data management classes and functions. Read in a deck of postcard images in
FITS files, match them up and return an array of tractor image data
structures. [Write out model images, residuals etc with sensible filenames.
Record key statistics in easily-parsed log files.]

'''

import numpy as np
import os,glob,string

import astrometry
from tractor import *

vb = 1

# ============================================================================

class deck:

    def __init__(self, folder):

        self.name = 'Deck of postcard images'
        self.dir = folder

        # Parse filenames for sci and wht images:
        self.inventory()

        # Read in data and parcel up as tractor.Images:
        self.riffle()

        return None

    def __str__(self):
        return '%s (%d in total)' % (self.name, self.number)

# ----------------------------------------------------------------------------
# Parse filenames for sci and wht images:

    def inventory(self):

        filenames = glob.glob(os.path.join(self.dir,'*.fits'))

        if vb: print "deck.inventory: found",len(filenames),"files: ",filenames

        # Break down file names. Naming convention: fruit_flavor.fits
        fruits = []
        flavors = []
        for filename in set(filenames):
           pieces = string.split(filename,'_')
           fruits.append(string.join(pieces[0:-1],'_'))
           flavors.append(string.split(pieces[-1],'.')[0])

        if len(set(flavors)) != 2:
           raise "deck.inventory: ERROR: expecting 2 flavors of datafile, got 0 or 1"

        if 'sci' not in set(flavors):
           raise "deck.inventory: ERROR: expecting some files to be xxx_sci.fits"

        for x in (set(flavors) - set(['sci'])):
           self.whttype = x

        self.number = len(set(fruits))
        self.scifile = []
        self.whtfile = []
        for fruit in set(fruits):
           x = fruit+'_sci.fits'
           if os.path.exists(x):
              self.scifile.append(x)
           else:
              self.scifile.append(None)
           x = fruit+'_'+self.whttype+'.fits'
           if os.path.exists(x):
              self.whtfile.append(x)
           else:
              self.whtfile.append(None)

        if vb:
           print "deck.inventory: packaged files into",self.number,"pair(s)"
           print "deck.inventory: 2 flavors of file found: sci and",self.whttype
           for i in range(self.number):
              print "deck.inventory:",i+1,"th pair:",[self.scifile[i],self.whtfile[i]]

        return

# ----------------------------------------------------------------------------
# Read in data and parcel up as tractor.Image instances:

    def riffle(self):

        return

# Hmm - need to build astrometry.net I think, to get the wcs
# functionality? Or make a ParamList wcs instead? This is better...
#     wcs = astrometry.util.util.Tan()
#     wcs.crval[0] = ra
#     wcs.crval[1] = dec
#     wcs.crpix[0] = W/2.
#     wcs.crpix[1] = H/2.
#     scale = width / float(W)
#     wcs.cd[0] = -scale
#     wcs.cd[1] = 0
#     wcs.cd[2] = 0
#     wcs.cd[3] = -scale
#     wcs.imagew = W
#     wcs.imageh = H
#
#     wcs = FitsWcs(wcs)
# ============================================================================

def PS1WCS(hdr):
      '''
      Return a WCS object initialised from a PS1 file header.
      '''
      t = astrometry.util.util.Tan()
      t.set_crpix(hdr['CRPIX1'], hdr['CRPIX2'])
      t.set_crval(hdr['CRVAL1'], hdr['CRVAL2'])
      t.set_cd(hdr['CDELT1'], 0., 0., hdr['CDELT2'])
      t.set_imagesize(hdr['NAXIS1'], hdr['NAXIS2'])
      
      return FitsWcs(t)
         
# ============================================================================

class LocalWCS(WCS):
      '''
      The "local" WCS -- useful when you need to work on the sky, but in
      small offsets from RA,Dec in arcsec. Initialisation is with the
      coordinates of the central pixel, which is set to be the origin of
      the "pixel coordinate" system.
      '''
      def __init__(self, pos, pixscale=1.0):
            self.x0 = 0.0
            self.y0 = 0.0
            self.ra = pos.ra
            self.dec = pos.dec
            self.pixscale = pixscale

      def cdAtPixel(self, x, y):
            return np.array([[1.0,0.0],[0.0,1.0]]) * self.pixscale / 3600.

      def hashkey(self):
            return ('LocalWCS', self.x0, self.y0, self.wcs)

      def __str__(self):
            return ('LocalWCS: x0,y0 %.3f,%.3f, WCS ' % (self.x0,self.y0)
                        + str(self.wcs))

      def setX0Y0(self, x0, y0):
            self.x0 = x0
            self.y0 = y0

      def positionToPixel(self, src, pos):
            # ok,x,y = self.wcs.radec2pixelxy(pos.ra, pos.dec)
            X = self.wcs.radec2pixelxy(pos.ra, pos.dec)
            if len(X) == 3:
                  ok,x,y = X
            else:
                  assert(len(X) == 2)
                  x,y = X
            return x-self.x0, y-self.y0

      def pixelToPosition(self, src, xy):
            (x,y) = xy
            r,d = self.wcs.pixelxy2radec(x + self.x0, y + self.y0)
            return RaDecPos(r,d)


# ============================================================================

if __name__ == '__main__':

# Basic test on lensfinder examples dir:

    folder = os.environ['LENSFINDER_DIR']+'/examples'
    postcards = deck(folder)

# ============================================================================
