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

class deck:

   def __init__(self, folder):

      self.name = 'Deck of postcard images'
      self.dir = folder

      # Parse filenames for sci and wht images:
      self.inventory()

#       # Read in data and parcel up as tractor.Images:
#       self.riffle()

      return None

# ----------------------------------------------------------------------------

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
# # Read in data and parcel up as tractor.Image instances:
# 
#    def riffle(self):
# 
#       # Loop over observations:
#       for i in range(self.number):
# 
# 
#       return
# 
# # ----------------------------------------------------------------------------
# 
#    def read_in_observation(self,i)
#       '''
#       Read in the ith sci and wht images, and make sure they have the right 
#       pixel values.
#       '''
#    
#    hdulist = pyfits.open(scifile)
#    sci = hdulist[0].data
#    hdr = hdulist[0].header
#    hdulist.close()
#    NX,NY = sci.shape
# 
#    hdulist = pyfits.open(varfile)
#    var = hdulist[0].data
#    hdulist.close()
#    assert sci.shape == var.shape
#    
#    mask = numpy.ones([NX,NY],dtype=np.int16)
#    mask[numpy.where(var == 0)] = 0
#    
#    # Convert var to wht, and find median uncertainty as well:
#    invvar = 1.0/var
#    # Assign zero weight to var=nan, var<=0:
#    invvar[var != var] = 0.0
#    invvar[var <= 0] = 0.0
#    
#    good = np.where(invvar > 0)
#    bad = np.where(invvar == 0)
#    
#    # Zero out sci image where wht is 0.0:
#    sci[bad] = 0.0
# 
#    assert(all(np.isfinite(sci.ravel())))
#    assert(all(np.isfinite(invvar.ravel())))
# 
#    # Rough estimates of background level and rms:
#    sciback = np.sqrt(np.median(var[good]))
#    scirms = np.sqrt(np.median(var[good]))
#    
#    # Report on progress so far:
#    if opt.verbose:
#       print 'Sci header:', hdr
#       print 'Read in sci image:', sci.shape, sci
#       print 'Read in var image:', var.shape, var
#       print 'Made mask image:', mask.shape, mask
#       print 'Useful variance range:', var[good].min(), var[good].max()
#       print 'Useful image median level:', sciback
#       print 'Useful image median pixel uncertainty:', scirms
# 
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

if __name__ == '__main__':

  
   if True:
   # Basic test on lenstractor examples dir:

      folder = os.environ['lenstractor_DIR']+'/examples'
      postcards = deck(folder)


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
