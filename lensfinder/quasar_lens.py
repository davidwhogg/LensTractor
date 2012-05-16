'''
This file is part of the LensFinder project.
Copyright 2012 David W. Hogg (NYU) and Phil Marshall (Oxford).

Description
-----------

Quasar lens model for lensed quasar detection. Given deflector parameters
(position,  Einstein radius, gammax, phix), and a source position, solve for
the 2,3 or 4 images, and make an image patch corresponding to this model.

To-do
-----
- debug, test

'''

from tractor.sdss_galaxy import DevGalaxy
from tractor import *
from gravitational_len import GravitationalLens

# ============================================================================
#
class EinsteinRadius(ScalarParam):
      def getName(self):
            return 'Einstein radius in pixels'

class ExternalShear(ParamList):
      def getName(self):
            return 'External shear, magnitude and angle'
      def getNamedParams(self):
            # gamma: shear magnitude, dimless
            # phi: deg, "E of N", 0=direction of increasing Dec, 90=direction of increasing RA
            return [('gamma', 0), ('phi', 1)]
      def hashkey(self):
            return ('ExternalShear',) + tuple(self.vals)
      def __repr__(self):
            return 'gamma=%g, phi=%g' % (self.gamma, self.phi)
      def __str__(self):
            return 'gamma=%.1f, phi=%.1f' % (self.gamma, self.phi)
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

class LensGalaxy(DevGalaxy):
      '''
      A LensGalaxy has mass, and emits light. Initialise with a position,
      brightness and shape for the Galaxy, and an Einstein radius and external
      shear for the mass distribution (assumed to be SIS+shear).
      '''

      def __init__(self, pos, brightness, shape, Rein, xshear):
            MultiParams.__init__(self, pos, brightness, shape, Rein, xshear)
            self.name = self.getName()

      def getName(self):
            return 'LensGalaxy'

      def getNamedParams(self):
            return [('pos', 0), ('brightness', 1), ('shape', 2),
                        ('Rein', 3), ('xshear', 4),]

      def getLensedImages(self,source)
            '''
            Unpack the parameters and pass to an instance of a
            gravitational lens, and ask for the image positions and fluxes for
            the given source. The source is a `PointSource` object.
            '''
            # Unpack the lens:
            lensRein = self.Rein
            lensgamma = self.gamma
            lensphi = np.deg2rad(self.phi) # lens solver expects radians.
            
            # Define a "trivial" coordinate system, centred on the lens, that
            # has 1 arcsec "pixels":
            tWCS = trivialWCS(self.pos) # Needs a new tractor class? Or inherit NullWCS
            lensposition = (0.0,0.0) # in trivial WCS
            
            # Unpack the source and convert position into trivial 
            # tangent-plane coordinates:
            sourceposition = tWCS.positionToPixel(source.getPosition())
            sourceflux = source.getBrightness()
            
            # Instantiate the gravitational lens:
            SISX = GravitationalLens(lensposition,lensRein,lensgamma,lensphi)
            
            # Solve for image positions and fluxes:
            fail,keepgoing = False,True
            imagepositions,fail = SISX.image_positions(sourceposition,keepgoing)
            imagefluxes = sourceflux*SISX.magnifications(imagepositions)

            # What to do when solver fails? 0.07% of the time...
            # Answer - deal with it. One or two images will be incorrectly 
            # merged or something, fine.
            if fail: pass 

            # Convert image positions back to sky:
            imagepositions = tWCS.pixelToPosition(imagepositions)
            
            return imagepositions, imagefluxes

# ============================================================================

class TrivialWCS(WCS):
	'''
	The "trivial" WCS -- useful when you need to work on the sky but in
      small offsets from RA,Dec.
	'''
	def __init__(self, pos, pixscale=1.):
            self.x0 = 0.0
            self.y0 = 0.0
            self.ra = pos.ra
            self.dec = pos.dec
		self.pixscale = pixscale
            # Set CD matrix?!
            

      # Somehow need to set wcs.cd... Is this how it is set?
	def cdAtPixel(self, x, y):
            cosdec = np.cos(self.dec*np.pi/180.0)
		return np.array([[cosdec,0.0],[0.0,1.0]]) * self.pixscale / 3600.
# 	def cdAtPixel(self, x, y):
# 		cd = self.wcs.get_cd()
# 		return np.array([[cd[0], cd[1]], [cd[2],cd[3]]])

	def hashkey(self):
		return ('TrivialWCS', self.x0, self.y0, self.wcs)

	def __str__(self):
		return ('TrivialWCS: x0,y0 %.3f,%.3f, WCS ' % (self.x0,self.y0)
				+ str(self.wcs))

	def setX0Y0(self, x0, y0):
		self.x0 = x0
		self.y0 = y0

	def positionToPixel(self, src, pos):
		# ok,x,y = self.wcs.radec2pixelxy(pos.ra, pos.dec)
		X = self.wcs.radec2pixelxy(pos.ra, pos.dec)
		if len(X) == 3:
			ok,x,y = X
		else:
			assert(len(X) == 2)
			x,y = X
		return x-self.x0, y-self.y0

	def pixelToPosition(self, src, xy):
		(x,y) = xy
		r,d = self.wcs.pixelxy2radec(x + self.x0, y + self.y0)
		return RaDecPos(r,d)


# ============================================================================

class QuasarLens(MultiParams):
       '''
       QuasarLens is a composite object consisting of a LensGalaxy [that has
       both light (Galaxy) and mass (GravitationalLens)], and a virtual Point
       Source behind it. The lensed image positions are determined by the
       lensing deflection; the magnified fluxes need to be perturbed in order
       to fit the data. Initialise with a LensGalaxy object, a PointSource
       object, and internally with four magnification perturbations (of which
       only 2, 3 or 4 are used in any given image configuration and which are
       initially set  to unity.
       '''

       def __init__(self, lensgalaxy, source):
            dmag = ParamList(np.zeros(4))
            # The next line fails if MultiParams.__init__ *copies*
            # self.lensgalaxy, self.source and self.dmag rather than points to
            # them...
            MultiParams.__init__(self, lensgalaxy, source, dmag)

       def getName(self):
               return 'QuasarLens'

       def getNamedParams(self):
               return [('lensgalaxy', 0), ('source', 1), ('dmag', 2)]

       def getModelPatch(self,img):
               '''
               Render the image of the quasarlens on the image grid provided.
               '''
               # Lens galaxy:
               patch = self.lensgalaxy.getModelPatch(img)
               # 2, 3 or 4 lensed images, from solving the lens equation.
               # Note: images are returned time-ordered:
               imagepositions, imagefluxes = self.lensgalaxy.getLensedImages(self.source)
               for i,(imageposition,imageflux) in enumerate(zip(imagepositions,imagefluxes)):
                  thisimageflux = imageflux * np.exp(dmag[i])
                  patch += PointSource(imageposition,thisimageflux).getModelPatch(img)

               return patch

       def getParamDerivatives(self, img, brightnessonly=False):
               pass


# ============================================================================

if __name__ == '__main__':

   pass
