'''
This file is part of the LensFinder project.
Copyright 2011 2012 2013 2014 David W. Hogg (NYU), Phil Marshall (KIPAC) and Adriano Agnello (UCLA).

Description
-----------

Nebula is a flexible model for a generic extended object, chosen to be the sum of np
point source components and ne extended source components. These extended sources - "NebulousGalaxies" - are just Dev galaxies with additional priors to avoid craziness.
'''

import numpy as np

import tractor
from tractor import galaxy

import lenstractor
from lenstractor import constants as const

# ============================================================================

class NebulousGalaxy(galaxy.DevGalaxy):

    '''
    A NebulousGalaxy is a Dev galaxy, initialized with a position, brightness,
    and shape, and given sensible Gaussian priors on all parameters.
    '''

    def __init__(self, pos, brightness, shape):
        tractor.MultiParams.__init__(self, pos, brightness, shape)
        self.name = self.getName()

    def __str__(self):
        return (self.getName() + ' at ' + str(self.pos) + ' with '
                               + str(self.pos.gpriors) + ', '
                               + ' with ' + str(self.brightness) + ' and '
                               + str(self.brightness.gpriors) + ', '
                               + ' and ' + str(self.shape) + ' with '
                               + str(self.shape.gpriors)
               )

    def getName(self):
        return 'NebulousGalaxy'

    def getNamedParams(self):
        return dict(pos=0, brightness=1, shape=2)

    def setPriors(self):
        # MAGIC: prior settings!
        self.pos.addGaussianPrior('ra',self.pos.ra,0.5*const.arcsec2deg*np.cos(self.pos.dec*const.deg2rad))
        self.pos.addGaussianPrior('dec',self.pos.dec,0.5*const.arcsec2deg)
        for band in self.brightness.order:
            self.brightness.addGaussianPrior(band,20.0,2.0)
        self.shape.addGaussianPrior('ee1', 0., 0.25)
        self.shape.addGaussianPrior('ee2', 0., 0.25)
        self.shape.addGaussianPrior('logre',np.log(0.5),0.3)
        return

# ============================================================================

if __name__ == '__main__':

    print "No tests defined, sorry."

# ============================================================================
