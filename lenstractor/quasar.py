'''
This file is part of the LensTractor project.
Copyright 2011 2012 2013 2014 David W. Hogg (NYU), Phil Marshall (KIPAC) and Adriano Agnello (UCLA).

Description
-----------

Quasars and Stars are models for point sources that are constrained by some priors to look like quasars and stars.
'''

import numpy as np

import tractor

import lenstractor
from lenstractor import constants as const

# ============================================================================

class Quasar(tractor.PointSource):

    '''
    A Quasar is a PointSource, initialized with a position and brightness,
    and given sensible Gaussian priors on all parameters.
    '''

    def __init__(self, pos, brightness):
        tractor.MultiParams.__init__(self, pos, brightness)
        self.name = self.getName()

    def __str__(self):
        return (self.getName() + ' at ' + str(self.pos) + ' with '
                               + str(self.pos.gpriors) + ', '
                               + ' with ' + str(self.brightness) + ' and '
                               + str(self.brightness.gpriors)
               )

    def getName(self):
        return 'Quasar'

    def getNamedParams(self):
        return dict(pos=0, brightness=1)

    def setPriors(self):
        # MAGIC: prior settings!
        self.pos.addGaussianPrior('ra',self.pos.ra,0.2*const.arcsec2deg*np.cos(self.pos.dec*const.deg2rad))
        self.pos.addGaussianPrior('dec',self.pos.dec,0.2*const.arcsec2deg)
        # Here is where we would insist on Quasar-like colors.
        for band in self.brightness.order:
            self.brightness.addGaussianPrior(band,20.0,2.0)
        return

# ============================================================================

class Star(tractor.PointSource):

    '''
    A Star is a PointSource, initialized with a position and brightness,
    and given sensible Gaussian priors on all parameters.
    '''

    def __init__(self, pos, brightness):
        tractor.MultiParams.__init__(self, pos, brightness)
        self.name = self.getName()

    def __str__(self):
        return (self.getName() + ' at ' + str(self.pos) + ' with '
                               + str(self.pos.gpriors) + ', '
                               + ' with ' + str(self.brightness) + ' and '
                               + str(self.brightness.gpriors)
               )

    def getName(self):
        return 'Star'

    def getNamedParams(self):
        return dict(pos=0, brightness=1)

    def setPriors(self):
        # MAGIC: prior settings!
        self.pos.addGaussianPrior('ra',self.pos.ra,2.0*const.arcsec2deg*np.cos(self.pos.dec*const.deg2rad))
        self.pos.addGaussianPrior('dec',self.pos.dec,2.0*const.arcsec2deg)
        # Here is where we would insist on Star-like colors.
        for band in self.brightness.order:
            self.brightness.addGaussianPrior(band,20.0,2.0)
        return

# ============================================================================

if __name__ == '__main__':

    print "No tests defined, sorry."

# ============================================================================
