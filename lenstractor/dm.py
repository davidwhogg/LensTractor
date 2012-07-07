# ============================================================================
'''
This file is part of the lenstractor project.
Copyright 2012 David W. Hogg (NYU) and Phil Marshall (Oxford).

Description
-----------

General-purpose data management classes and functions: 
* Order a pile of FITS files into scifiles and matching varfiles
* Read in a deck of postcard images in FITS files and return an 
  array of tractor image data structures. 
'''

import numpy as np
import os,glob,string,pyfits

from astrometry.util import util
import tractor
import lenstractor

# ============================================================================
# Parse filenames for sci and wht images:

def Riffle(filenames,vb=False):

   if vb: print "Looking at",len(filenames),"files: ",filenames

   # Break down file names. Naming convention: fruit_flavor.fits
   fruits = []
   flavors = []
   for filename in set(filenames):
      pieces = string.split(filename,'_')
      fruits.append(string.join(pieces[0:-1],'_'))
      flavors.append(string.split(pieces[-1],'.')[0])

   if len(set(flavors)) != 2:
      raise "ERROR: expecting 2 flavors of datafile, got 0 or 1"

   if 'sci' not in set(flavors):
      raise "ERROR: expecting at least some files to be xxx_sci.fits"

   for x in (set(flavors) - set(['sci'])):
      whttype = x

   number = len(set(fruits))
   scifiles = []
   whtfiles = []
   for fruit in set(fruits):
      x = fruit+'_sci.fits'
      if os.path.exists(x):
         scifiles.append(x)
      else:
         scifiles.append(None)
      x = fruit+'_'+whttype+'.fits'
      if os.path.exists(x):
         whtfiles.append(x)
      else:
         whtfiles.append(None)

   if vb:
      print "Riffled files into",number,"pair(s)"
      print "2 flavors of file found: sci and",whttype
      for i in range(number):
         print "   ",i+1,"th pair:",[scifiles[i],whtfiles[i]]

   return scifiles,whtfiles

# ============================================================================
# Read in data ad organise into Tractor Image objects.
# Some of this is survey specific: subroutines to be stored in $survey.py.

def Deal(scifiles,varfiles,SURVEY='PS1',vb=False):

   images = []
   bands = []
   epochs = []
   
   for scifile,varfile in zip(scifiles,varfiles):
      
      name = scifile.replace('_sci.fits','')
      if vb: 
         print " "
         print "Making Tractor image from "+name+"_*.fits:"

      # Read in sci and wht images. Note assumptions about file format:
      sci,invvar,hdr = Read_in_data(scifile,varfile,vb)
      
      # Initialize a PSF object (single Gaussian by default), first 
      # getting FWHM from somewhere. Start with FWHM a little small, 
      # then refine it:

      if SURVEY=='PS1':
         FWHM = lenstractor.PS1_IQ(hdr)
      else:
         Raise("Unrecognised survey %s" % SURVEY)
      if vb: print "  PSF FWHM =",FWHM,"pixels"

      # MAGIC 0.8 shrinkage factor:
      psf = Initial_PSF(0.8*FWHM)
      if vb: print psf

      # Now get the photometric calibration from the image header.

      if SURVEY=='PS1':
         band,photocal = lenstractor.PS1_photocal(hdr)
      else:
         Raise("Unrecognised survey %s" % SURVEY)
      if vb: print photocal
      bands.append(band)
      epochs.append(lenstractor.PS1_epoch(hdr))
      
      # Set up sky to be varied:
      sky = tractor.ConstantSky(0.0)
      if vb: print sky

      # Get WCS from FITS header:
      wcs = lenstractor.PS1WCS(hdr)
      if vb: print wcs

      # Make a tractor Image object out of all this stuff, and add it to the array:
      images.append(tractor.Image(data=sci, invvar=invvar, name=name,
				    psf=psf, wcs=wcs, sky=sky, photocal=photocal))

   # Figure out the unique band names and epochs:
   uniqbands = np.unique(np.array(bands))
   if vb: 
      print " "
      print "Read in",len(images),"image datasets"
      print "  in",len(uniqbands),"bands:",uniqbands
      print "  at",len(epochs),"epochs"
      print " "

   return images,np.array(uniqbands)
      
# ============================================================================
# Read in sci and wht images. Note assumptions about file format:

def Read_in_data(scifile,varfile,vb=False):

   hdulist = pyfits.open(scifile)
   sci = hdulist[0].data
   hdr = hdulist[0].header
   hdulist.close()
   NX,NY = sci.shape

   hdulist = pyfits.open(varfile)
   var = hdulist[0].data
   hdulist.close()
   assert sci.shape == var.shape

   # Convert var to wht, and find median uncertainty as well:
   invvar = 1.0/var
   # Assign zero weight to var=nan, var<=0:
   invvar[var != var] = 0.0
   invvar[var <= 0] = 0.0

   good = np.where(invvar > 0)
   bad = np.where(invvar == 0)

   # Zero out sci image where wht is 0.0:
   sci[bad] = 0.0

   assert(all(np.isfinite(sci.ravel())))
   assert(all(np.isfinite(invvar.ravel())))

   # Very rough estimates of background level and rms, never used:
   sciback = np.median(sci[good])
   scirms = np.sqrt(np.median(var[good]))

   # Report on progress so far:
   if vb:
      # print 'Sci header:', hdr
      print 'Read in sci image:', sci.shape #, sci
      print 'Read in var image:', var.shape #, var
      # print 'Made mask image:', mask.shape, mask
      print 'Useful variance range:', var[good].min(), var[good].max()
      print 'Useful image median level:', sciback
      print 'Useful image median pixel uncertainty:', scirms
      
   return sci,invvar,hdr
      
# ============================================================================
# Initialize a PSF object - by default, a single circularly symmetric Gaussian 
# defined on same grid as sci image:

def Initial_PSF(FWHM,double=False):

   # NB. FWHM of PSF is given in pixels.

   if not double:
      # Single Gaussian default:
      w = np.array([1.0])                      # amplitude at peak
      mu = np.array([[0.0,0.0]])               # centroid position in pixels 
      var = (FWHM/2.35)**2.0
      cov = np.array([[[var,0.0],[0.0,var]]])  # pixels^2, covariance matrix

   else:
      # Double Gaussian alternative:
      w = np.array([1.0,1.0])                                        
      mu = np.array([[0.0,0.0],[0.0,0.0]])                           
      var = (FWHM/2.35)**2.0
      cov = np.array([[[1.0,0.0],[0.0,1.0]],[[var,0.0],[0.0,var]]])  
   
   return tractor.GaussianMixturePSF(w,mu,cov)

# ============================================================================
if __name__ == '__main__':

   if True:
   
   # Basic test on lenstractor examples dir:
    
      folder = os.environ['LENSTRACTOR_DIR']+'/examples'
      inputfiles = glob.glob(os.path.join(folder,'*.fits'))

      scifiles,varfiles = riffle(inputfiles)

      

