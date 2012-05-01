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

# import tractor
# BUG: this fails because I don't have astrometry.sdss
#      Very likely there are more problems like this!

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

# ============================================================================

if __name__ == '__main__':

# Basic test on lensfinder examples dir:

    folder = os.environ['LENSFINDER_DIR']+'/examples'
    postcards = deck(folder)
   
# ============================================================================
