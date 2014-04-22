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
import os,glob,string,pyfits,subprocess

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

   if len(set(flavors)) > 2:
       raise ValueError("ERROR: expecting 1 or 2 flavors of datafile, got more")
   elif len(set(flavors)) == 0:
       raise ValueError("ERROR: expecting 1 or 2 flavors of datafile, got none")

   if 'sci' not in set(flavors):
       raise ValueError("ERROR: expecting at least some files to be xxx_sci.fits")

   if len(set(flavors)) == 1:
       whttype = 'No-one_will_ever_choose_this_flavor'   
   else:
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
      if len(set(flavors)) == 1:     
          print "Only 1 flavor of file found, sci"
      else:
          print "2 flavors of file found: sci and",whttype
      for i in range(number):
         print "   ",i+1,"th pair:",[scifiles[i],whtfiles[i]]

   return scifiles,whtfiles

# ============================================================================
# Read in data and organise into Tractor Image objects.
# Some of this is survey specific: subroutines to be stored in $survey.py.

def Deal(scifiles,varfiles,SURVEY='PS1',vb=False):

   images = []
   bands = []
   epochs = []
   total_mags = []
   
   for scifile,varfile in zip(scifiles,varfiles):
      
      name = scifile.replace('_sci.fits','')
      if vb: 
         print " "
         print "Making Tractor image from "+name+"_*.fits:"

      # Read in sci and wht images. Note assumptions about file format:
      sci,invvar,hdr,total_flux = Read_in_data(scifile,varfile,vb)
      
      if total_flux == 0.0:
         print "No flux found in image from "+scifile
         print "Skipping to next image!"
         continue
      
      # Initialize a PSF object (single Gaussian by default), first 
      # getting FWHM from somewhere. Start with FWHM a little small, 
      # then refine it:

      if SURVEY=='PS1':
         try:
             FWHM = lenstractor.PS1_IQ(hdr)
         except:
             FWHM = 1.4
      elif SURVEY=='KIDS':
         FWHM = lenstractor.KIDS_IQ(hdr)
      elif SURVEY=='SDSS':
         FWHM = lenstractor.SDSS_IQ(hdr)
      else:
         raise ValueError('Unrecognised survey name '+SURVEY)
      if vb: print "  PSF FWHM =",FWHM,"pixels"

      # MAGIC shrinkage factor:
      shrink = 0.8
      psf = Initial_PSF(shrink*FWHM)
      if vb: print psf

      # Now get the photometric calibration from the image header.

      if SURVEY=='PS1':
         try:
             band,photocal = lenstractor.PS1_photocal(hdr)
         except:
             band,photocal = lenstractor.SDSS_photocal(hdr)
      elif SURVEY=='KIDS':
         band,photocal = lenstractor.KIDS_photocal(hdr)
      elif SURVEY=='SDSS':
         band,photocal = lenstractor.SDSS_photocal(hdr)
      else:
         print "Unrecognised survey name "+SURVEY+", assuming SDSS"
         band,photocal = lenstractor.SDSS_photocal(hdr)

      if vb: print photocal
      bands.append(band)
      if SURVEY=='PS1':
         try:
             epochs.append(lenstractor.PS1_epoch(hdr))
         except:
             epochs.append(lenstractor.SDSS_epoch(hdr))
      elif SURVEY=='KIDS':
         epochs.append(lenstractor.KIDS_epoch(hdr))
      elif SURVEY=='SDSS':
         epochs.append(lenstractor.SDSS_epoch(hdr))
      
      # Use photocal to return a total magnitude:
      total_mag = photocal.countsToMag(total_flux)
      if vb: print "Total brightness of image (mag):",total_mag
      total_mags.append(total_mag)

      # Set up sky to be varied:
      median = np.median(sci[invvar > 0])
      sky = tractor.ConstantSky(median)
      delta = 0.1*np.sqrt(1.0/np.sum(invvar))
      assert delta > 0
      sky.stepsize = delta
      if vb: print sky

      # Get WCS from FITS header:
      if SURVEY=='PS1':
         try:
             wcs = lenstractor.PS1WCS(hdr)
         except:
             wcs = lenstractor.SDSSWCS(hdr)
      elif SURVEY=='KIDS':
         wcs = lenstractor.KIDSWCS(hdr)
      else:
         try:
             wcs = lenstractor.SDSSWCS(hdr)
         except:
             wcs = lenstractor.SDSSWCS(hdr)
      # if vb:
      #    print wcs

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

   return images,np.array(total_mags),np.array(bands)
      
# ============================================================================
# Read in sci and wht images. Note assumptions about file format:

def Read_in_data(scifile,varfile,vb=False):

   hdulist = pyfits.open(scifile)
   sci = hdulist[0].data
   hdr = hdulist[0].header
   hdulist.close()
   NX,NY = sci.shape

   if (varfile is not None): 
       hdulist = pyfits.open(varfile)
       var = hdulist[0].data
       hdulist.close()
   else:
   # Make a var image from the sci image...
       background = np.median(sci)
       diffimage = sci - background
       variance = np.median(diffimage*diffimage)
       var = diffimage + variance
       # Ensure positivity:
       var[var <= 0] = variance 

   # Check image sizes...
   assert sci.shape == var.shape

   # Convert var to wht, and find median uncertainty as well:
   invvar = 1.0/var
   # Assign zero weight to var=nan, var<=0:
   invvar[var != var] = 0.0
   invvar[var <= 0] = 0.0

   bad = np.where(invvar == 0)
   # Zero out sci image where wht is 0.0:
   sci[bad] = 0.0

   assert(all(np.isfinite(sci.ravel())))
   assert(all(np.isfinite(invvar.ravel())))

   # Measure total flux in sci image:
   total_flux = np.sum(sci)

   # Report on progress so far:
   if vb:
      print 'Science image:', sci.shape #, sci
      print 'Total flux:', total_flux
      print 'Variance image:', var.shape #, var

   if total_flux != 0.0:
      # Very rough estimates of background level and rms, never used:
      good = np.where(invvar > 0)
      sciback = np.median(sci[good])
      scirms = np.sqrt(np.median(var[good]))
      if vb:
         print 'Useful variance range:', var[good].min(), var[good].max()
         print 'Useful image median level:', sciback
         print 'Useful image median pixel uncertainty:', scirms
      
   return sci,invvar,hdr,total_flux
      
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

