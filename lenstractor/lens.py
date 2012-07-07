'''
This file is part of the lenstractor project.
Copyright 2012 David W. Hogg (NYU) and Phil Marshall (Oxford).

Description
-----------

Lens objects for lens detection. Given deflector parameters
(position,  Einstein radius, gammax, phix), and a source position, solve for
the 2,3 or 4 images of a point source lens, and make an image patch 
corresponding to this model.

To-do
-----
- debug, test

'''

import numpy as np

from astrometry.util import util

import tractor
import lenstractor


# ============================================================================
# Gravitational Lens parameters need to be "Params" objects of some kind.

class EinsteinRadius(tractor.ScalarParam):
      def getName(self):
            return 'Einstein radius'

class ExternalShear(tractor.ParamList):
      def getName(self):
            return 'External shear'
      def getNamedParams(self):
            # gamma: shear magnitude, dimless
            # phi: deg, "E of N", 0=direction of increasing Dec, 90=direction of increasing RA
            return dict(gamma=0, phi=1)
      def hashkey(self):
            return ('ExternalShear',) + tuple(self.vals)
      def __repr__(self):
            return 'gamma=%g, phi=%g' % (self.gamma, self.phi)
      def __str__(self):
            return '%s: gamma=%.1f, phi=%.1f' % (self.getName(), self.gamma, self.phi)
      def copy(self):
            return ExternalShear(*self.vals)
      def getParamNames(self):
            return ['gamma', 'phi']
      def getStepSizes(self, img):
            return [ 0.1, 1.0 ]

# ============================================================================
# Composite object, consisting of:
#  - deV galaxy
#  - SIS + shear mass distribution to act as a gravitational lens

class LensGalaxy(tractor.sdss_galaxy.DevGalaxy):
      '''
      A LensGalaxy has mass, and emits light. Initialise with a position,
      brightness and shape for the Galaxy (which is assumed to have a De 
      Vaucoleurs surface brightness profile), and an Einstein radius and 
      external shear for the mass distribution (assumed to be SIS).
      '''

      def __init__(self, pos, brightness, shape, Rein, xshear):
            tractor.MultiParams.__init__(self, pos, brightness, shape, Rein, xshear)
            self.name = self.getName()

      def __str__(self):
            return (self.getName() + ' at ' + str(self.pos)
                                   + ' with ' + str(self.brightness) 
                                   + ' and ' + str(self.shape)
                                   + ' and ' + str(self.Rein)
                                   + ' and ' + str(self.xshear))
      
      def getName(self):
            return 'LensGalaxy'

      def getNamedParams(self):
            return dict(pos=0, brightness=1, shape=2, Rein=3, xshear=4)

      def getLensedImages(self,source):
            '''
            Unpack the parameters and pass to an instance of a
            gravitational lens, and ask for the image positions and fluxes for
            the given source. The source is a `PointSource` object.
            '''
            # Unpack the lens:
            lensRein = self.Rein.getValue()
            lensgamma = self.xshear[0]
            lensphi = np.deg2rad(self.xshear[0]) # lens solver expects radians.
            
            # Define a "trivial" coordinate system, centred on the lens, that
            # has 1 arcsec "pixels":
            lenswcs = lenstractor.LensPlaneWCS(self.pos) # Trivial tangent plane wcs, 1" pixels, N up
            lenspixelpos = (0.0,0.0)                    # in Lens Plane WCS
            
            # Unpack the source and convert position into trivial 
            # tangent-plane coordinates:
            sourcepixelpos = np.array(lenswcs.positionToPixel(source.getPosition()))
            
            # Instantiate the gravitational lens:
            SISX = lenstractor.GravitationalLens(lenspixelpos,lensRein,lensgamma,lensphi)
            
            # Solve for image positions and fluxes:
            fail,keepgoing = False,True
            imagepixelpos,fail = SISX.image_positions(sourcepixelpos,keepgoing)
            imagemagnifications = SISX.magnifications(imagepixelpos)

            # What to do when solver fails? 0.07% of the time...
            # Answer - deal with it. One or two images will be incorrectly 
            # merged or something, fine.
            if fail: pass 

            # Convert image positions back to sky:
            imagepositions = [lenswcs.pixelToPosition(p[0],p[1]) for p in imagepixelpos]
                        
            return imagepositions, imagemagnifications

# ============================================================================

class PointSourceLens(tractor.MultiParams):
       '''
       PointSourceLens is a composite object consisting of a LensGalaxy [that has
       both light (Galaxy) and mass (GravitationalLens)], and a virtual Point
       Source behind it. The lensed image positions are determined by the
       lensing deflection; the magnified fluxes need to be perturbed in order
       to fit the data. Initialise with a LensGalaxy object, a PointSource
       object, and internally with four magnification perturbations (of which
       only 2, 3 or 4 are used in any given image configuration and which are
       initially set  to unity.
       '''

       def __init__(self, lensgalaxy, pointsource):
#             dmag = ParamList(np.zeros(4))
            # The next line fails if MultiParams.__init__ *copies*
            # self.lensgalaxy, self.pointsource and self.dmag rather than 
            # points to them...
#             MultiParams.__init__(self, lensgalaxy, pointsource, dmag)
            tractor.MultiParams.__init__(self, lensgalaxy, pointsource)
            
            # Create 4 local cached PointSource instances for the purpose of 
            # patch-making later:
            self.pointsourcecache = [pointsource.copy() for i in range(4)]
          
            return

       def __str__(self):
            return (self.getName() + ' comprising a ' + str(self.lensgalaxy)
                                   + ' and a ' + str(self.pointsource))
                                   
       def getName(self):
               return 'PointSourceLens'

       def getNamedParams(self):
#                return dict(lensgalaxy=0, pointsource=1, dmag=2)
               return dict(lensgalaxy=0, pointsource=1)

       def getModelPatch(self,img):
               '''
               Render the image of the PointSourceLens on the image grid provided.
               '''
               # Lens galaxy:
               patch = self.lensgalaxy.getModelPatch(img)
               
               # Solve the lens equation to get the image positions and fluxes.
               # Note: images are returned time-ordered:
               imagepositions, imagemagnifications = self.lensgalaxy.getLensedImages(self.pointsource)
                              
               # Add point image patches to the patch, applying dmags:
               for i,(imageposition,imagemagnification) in enumerate(zip(imagepositions,imagemagnifications)):
                  # Recall: pointsourcecache is a list of 4 pointsource instances, to be pointed at.
                  PS = self.pointsourcecache[i]
                  PS.setPosition(imageposition)
                  PS.setBrightness(self.pointsource.getBrightness()*imagemagnification)
                  # This is brittle - if the Patch is entirely outside the FoV, getModelPatch returns None, 
                  # which cannot be added...
                  patch += PS.getModelPatch(img)

               return patch

#        def getParamDerivatives(self, img, brightnessonly=False):
#                # Basic parameter derivatives by finite differencing:
#                pars0 = self.getParams()
#                patch0 = self.getModelPatch(img)
#                derivs = []
#                for i,(step,name) in enumerate(zip(self.getStepSizes(), self.getParamNames())):
#                        print 'Img', img.name, 'deriv', i, name
#                        oldval = self.setParam(i, pars0[i] + step)
#                        patchi = self.getModelPatch(img)
#                        self.setParam(i, oldval)
#                        dpatch = (patchi - patch0) / step
#                        derivs.append(Patch(0, 0, dpatch))
#                return derivs


# ============================================================================

def LensPlaneWCS(pos):
      '''
      The "local" WCS -- useful when you need to work on the sky, but in
      small offsets from RA, Dec in arcsec. Initialisation is with the
      coordinates of the central pixel, which is set to be the origin of
      the "pixel coordinate" system. Return a WCS object.
      '''
      
      onearcsec = 1.0/3600.0
      
      return tractor.FitsWcs(util.Tan(pos.ra,pos.dec,1.0,1.0,-onearcsec,0.0,0.0,onearcsec,0,0))


# ============================================================================

if __name__ == '__main__':

   pass
