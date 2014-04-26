'''
This file is part of the lenstractor project.
Copyright 2012 David W. Hogg (NYU) and Phil Marshall (Oxford).

Description
-----------

Model class, to package up srcs and enable initialization.

To-do
-----
- debug, test

'''

import numpy as np

import tractor
import lenstractor

# ============================================================================

class Model():
    '''
    PURPOSE
      Define a named model that can be initialized, passed to a tractor etc.

    COMMENTS

    CONTENTS
      srcs             A list of tractor.sources
      name             etc

    OUTPUTS

    BUGS

    HISTORY
      2014-04-17       Started Marshall & Agnello (UCSB)
    '''
# ----------------------------------------------------------------------------
    
    def __init__(self,name,srcs=None,vb=True):
    
        self.name = name
        self.flavor = name[0:6]
        if self.flavor == 'Nebula':
           # How many point sources?
           self.K = int(name[6:7])
        
        if srcs == None:
            self.srcs = []
        
        self.vb = vb
        
        return None
        
# ----------------------------------------------------------------------------
    
    def __str__(self):
        return '%s' % (self.name)

# ----------------------------------------------------------------------------
    
    def initialize(self,template,position=None,SED=None):
                
        if template == 'from_scratch':
                        
            assert position is not None
            assert SED is not None
            
            print "Initializing",self.name,"model from scratch..."
           
            if self.flavor == 'Nebula':
                self.create_Nebula(position,SED)
            else:
                self.create_Lens(position,SED)

        else:
            
            # Initialization from template:

            print "Initializing",self.name,"model from",template.name," template..."
            pass
        
        return None
            
# ----------------------------------------------------------------------------
    
    def create_Nebula(self,position,SED):
                
        # Start with an exponential galaxy at the object centroid:
        galpos = position
        # Make it fainter than the point sources, by a couple of magnitudes
        fudge = 2
        galSED = SED.copy() + 2.5*np.log10(5*fudge)
        # Some standard shape and size parameters:
        re = 0.5    # arcsec
        q = 0.8     # axis ratio
        theta = 0.0 # degrees
        galshape = tractor.sdss_galaxy.GalaxyShape(re,q,theta)
        # Package up:
        nebulousgalaxy = tractor.sdss_galaxy.ExpGalaxy(galpos,galSED,galshape)
        if self.vb: print nebulousgalaxy
        self.srcs.append(nebulousgalaxy)

        # Now add the point sources:
        for i in range(self.K):
            # Small random offsets from nebula centre:
            e = 0.2 # arcsec
            dx,dy = e*np.random.randn(2)/3600.0
            starpos = position.copy() + tractor.RaDecPos(dx,dy)
            # Package up:
            star = tractor.PointSource(starpos,SED.copy())
            if self.vb: print star
            self.srcs.append(star)
            
        return

# ----------------------------------------------------------------------------
    
    def create_Lens(self,position,SED):
                
        # Start with a source to be lensed:
        xs = position.copy()
        ms = SED.copy() + 2.5*np.log10(40.0)
        pointsource = tractor.PointSource(xs,ms)
        if self.vb: print pointsource

        # Figure out a suitable initial lens potential:
        thetaE = lenstractor.EinsteinRadius(0.2) # arcsec. Start small, Adri?
        gamma = 0.2  # to make a quad
        phi   = 45.0 # deg
        xshear = lenstractor.ExternalShear(gamma,phi)
        # Add the lens galaxy light:
        xd = position.copy()
        md = SED.copy() + 2.5*np.log10(10.0)
        re = 0.5  # arcsec
        q = 0.8   # axis ratio
        theta = -phi # Note how mass/light misalignment is enabled.
        galshape = tractor.sdss_galaxy.GalaxyShape(re,q,theta)
        lensgalaxy = lenstractor.LensGalaxy(xd,md,galshape,thetaE,xshear)
        if self.vb: print lensgalaxy

        self.srcs.append(lenstractor.PointSourceLens(lensgalaxy, pointsource))

        return
        
# ============================================================================

if __name__ == '__main__':

   pass
