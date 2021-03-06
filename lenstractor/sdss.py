# ============================================================================

'''
This file is part of the lenstractor project.
Copyright 2014 Phil Marshall (KIPAC).

Description
-----------

Data management classes and functions for the SDSS survey.

* Read in a deck of postcard images in FITS files, match them up and return an 
  array of tractor image data structures. 
* Manipulate SDSS FITS image headers, WCS and photo cal items.
'''

import numpy as np
import os,glob,string

from astrometry.util import util
import tractor
from tractor import sdss as st
import lenstractor

vb = True

# ============================================================================

def SDSSWCS(hdr):
      '''
      Return a WCS object initialised from a SDSS file header.
      '''
      t = util.Tan()
      t.set_crpix(hdr['CRPIX1'], hdr['CRPIX2'])
      t.set_crval(hdr['CRVAL1'], hdr['CRVAL2'])
      t.set_cd(hdr['CD1_1'], 0., 0., hdr['CD2_2'])
      t.set_imagesize(hdr['NAXIS1'], hdr['NAXIS2'])
      
      return tractor.FitsWcs(t)

# (Adri's version)
# def SDSSWCS(hdr):
#      '''
#      Return a WCS object initialised from a PS1 file header.
#      WARNING: PC matrix not being used, determinant may be wrong sign...
#      Need to check with literature image of H1413
#      '''
#      t = util.Tan()
#      t.set_crpix(hdr['CRPIX1'], hdr['CRPIX2'])
#      t.set_crval(hdr['CRVAL1'], hdr['CRVAL2'])
#      t.set_cd(hdr['CD1_2'], 0., 0., hdr['CD2_2'])
#      t.set_imagesize(hdr['NAXIS1'], hdr['NAXIS2'])
#      
#      return tractor.FitsWcs(t)
#

# ============================================================================

def SDSS_IQ(hdr):
   """
   Return the image quality (FWHM) in pixels given a FITS file header.
   """
   # TODO: Allow for some FITS images that do have PSF information?
   
   FWHM = 1.4 # arcsec
   # Need it in pixels:
   determinant = hdr['CD1_1']*hdr['CD2_2'] - hdr['CD1_2']*hdr['CD2_1']
   pixscale = np.sqrt(np.abs(determinant))
   print pixscale
   pixscale *= 3600.0
   print "pixel-scale = ",pixscale
   FWHM /= pixscale
    
   return FWHM


# ============================================================================

class SDSSMagsPhotoCal(tractor.BaseParams):
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
            return 'SDSSMagsPhotoCal'
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
            return 10.**(0.4 * (self.zpt - mag))

      def countsToMag(self, counts):
            # Returns cheap scalar magnitude
            if counts <= 0.0:
                  return 99.0
            # Assume zeropoint is the whole story:
            return self.zpt - 2.5*np.log10(counts)


def SDSS_photocal(hdr):
   """
   Return a PS1MagsPhotoCal object given an SDSS FITS file header.
   Ignore "NMGY" keyword, and copy code found in tractor/tractor/sdss.py 
   class SdssNanomaggiesPhotoCal. Resulting mags seem to be sensible...
   """

   band = hdr['FILTER'][0]
   zpt= 22.5
   
   photocal = lenstractor.SDSSMagsPhotoCal(zpt,band)

   return band,photocal

# (Adri's version)
# def SDSS_photocal(hdr):
#
#   band = hdr['FILTER'][0]
#   zpt = hdr['NMGY']
#   zpt= -2.5*np.log10(zpt*3.631*10**(-6))
#   photocal = lenstractor.SDSSMagsPhotoCal(zpt,band)
#   return band,photocal

# ============================================================================
# 
# def SDSS_photocal(hdr):
#    """
#    Return a SDSS PhotoCal object given an SDSS FITS file header.
#    """
# 
#    bandname = hdr['FILTER'][0]
#    photocal = tractor.SdssNanomaggiesPhotoCal(bandname)
# # AttributeError: 'module' object has no attribute 'SdssNanomaggiesPhotoCal'
#    return band,photocal
# 
# ============================================================================

def SDSS_epoch(hdr):
   """
   Return the observation date (not MJD!) given a FITS file header.
   """
   # TO DO: convert to MJD.
   
   return "%s" % hdr['DATE-OBS']

# ============================================================================
# Return scales to use in Plot_state() of plotting.py
# image is an Image instance from tractor.

def SDSS_imshow_settings(image, chi):
    # The image is inverted when plotting, therefore
    # 'min' and 'max' are inverted here as well. 
    theimage = image.getImage()
    ima = dict(interpolation='nearest', origin='lower',
               vmin=-theimage.max(), vmax=-theimage.min())

    # The values in the differential chi images are in the range
    # of 0.0001, while the plotting routines expects them to be
    # within the range of 1. The values are apparently supposed to
    # be 'sigmas' but this cannot be the case for SDSS.
    # The assumption that the chi images represent sigmas might
    # be made at other places in the code as well, which might
    # be a problem.
    chia = dict(interpolation='nearest', origin='lower',
                vmin=-chi.max(), vmax=-chi.min())

    psfa = dict(interpolation='nearest', origin='lower')
 
    return (ima, chia, psfa)

# ============================================================================
# Get SDSS data, dstndstn style!

def getSDSSdata(rcf,roi,datadir,vb=False):
    '''
    INPUTS
        rcf          A list of number, specifying the SDSS run, camcol and field ids, ra and dec
        roi          Width of region of interest (in arcsec)
        datadir      Name of directory where fits images will be written to (for posterity)
        vb           Verbose operation?

    OUTPUTS
        images       A list of tractor image objects, one per image
        centroids    List of flux centroids, one per image
        total_mags   List of total magnitudes, one per image
        bands        List of filter names, one per image
    '''
    
    run,camcol,field,ra,dec = rcf[0],rcf[1],rcf[2],rcf[3],rcf[4]
    roipix = int(roi/0.396)
    geometry = (ra, dec, 0.5*roipix)
    SDSSbands = ['u','g','r','i','z']
    bands = []
    images = []
    centroids = []
    total_mags = []
    
    if vb:
        print "Querying SDSS skyserver for data at ra,dec = ",ra,dec
        print "  in run, camcol, field = ",[run,camcol,field]
        print "Looking to make cutouts that are ",roi," arcsec (",roipix," pixels) across"
    
    # Ideally we would have first grabbed SDSS run, camcol and field id:
    #   rcf = radec_to_sdss_rcf(ra,dec,radius=math.hypot(fieldradius,13./2.),tablefn="dr9fields.fits")
    # but, we do not have this table. So need an rcf passed in as an argument.

    for band in SDSSbands:
    
        if vb:
            print " "
            print "--------"+band+"-band--------"
        
        # Get SDSS image:
        image,info = st.get_tractor_image_dr9(run, camcol, field, band, roiradecsize=geometry)

        bands.append(band)
        images.append(image)

        # Compute flux centroid for this image (dummy for now):
        centroid = tractor.RaDecPos(ra,dec)
        centroids.append(centroid)

        # Compute total magnitude for this image:
        background = np.median(image.data)
        diffimage = image.data - background
        total_flux = np.sum(diffimage)
        total_mag = image.photocal.countsToMag(total_flux)
        total_mags.append(total_mag)
        
        if vb: 
            print "Got pixel data in the "+band+" band"
            print "Image size: ",image.data.shape,"pixels"
            print "Total brighness of image (mags): ",total_mag


    return images,centroids,np.array(total_mags),np.array(bands)

# ============================================================================
