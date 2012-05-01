from tractor import sdss_galaxy as stgal
from tractor import *

# ============================================================================

class LensGalaxy(MultiParams):

       def __init__(self, light, mass):
               MultiParams.__init__(self, light, mass)
               self.name = self.getName()

       def getName(self):
               return 'Lens Galaxy'

       def getNamedParams(self):
               return [('light', 0), ('mass', 1)]

# ============================================================================

class LensedQuasar(MultiParams):

       def __init__(self, lens, source):
               MultiParams.__init__(self, mass, source)
               self.name = self.getName()

       def getName(self):
               return 'LensedQuasar'

       def getNamedParams(self):
               return [('mass', 0), ('source', 1)]

       def hashkey(self):
               # important to implement -- but MultiParams does it for you?
               pass

       def getModelPatch(self, img, px=None, py=None):
               # real action here -- return a Patch object
               # -this might be the only place you create the 4 PointSource 
               # images.
               pass

       def getParamDerivatives(self, img, brightnessonly=False):
               # pass


# ============================================================================

if __name__ == '__main__':

# Make a galaxy, with light and mass:

       light = stgal.CompositeGalaxy(pos, brexp, shexp, brdev, shdev)
       mass = ....
       lens = LensGalaxy(light, mass)

       source = PointSource(pos, bright)
       images = LensedQuasar(lens, source)

       lq.getParams()
       # --> list of ~20 param values (floating-point vals) - from where?

       # Have a look at CompositeGalaxy
