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
        return (self.getName() + ' at ' + str(self.pos)
                               + ' with ' + str(self.brightness)
                               + ' and ' + str(self.shape)
               )

    def getName(self):
            return 'NebulousGalaxy'

    def getNamedParams(self):
            return dict(pos=0, brightness=1, shape=2)

# ============================================================================

if __name__ == '__main__':

    print "No tests defined, sorry."

# ============================================================================
