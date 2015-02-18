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
import pylab as plt

import tractor
import lenstractor

deg2rad = 0.017453292
resetre = 0.5 # effective radius for the galaxy reset, 0.5 for PS1 and 1.0 for SQLS

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
#        else:
#           self.K = 5

        if srcs == None:
            self.srcs = []

        self.vb = vb

        return None

# ----------------------------------------------------------------------------

    def __str__(self):
        return '%s' % (self.name)

    def copy(self):
        copy = lenstractor.Model(self.name,vb=self.vb)
        for src in self.srcs:
            copy.srcs.append(src.copy())
        return copy

# ----------------------------------------------------------------------------

    def initialize(self,template,position=None,SED=None):

        if template == 'from_scratch':

            # Are we initializing positions from a catalog?
            if isinstance(position,basestring):

                if self.vb: print "Initializing",self.name,"model from catalog..."

                self.manual = True
                position = self.get_positions_from(position,SED)

            else:
                assert position is not None
                assert SED is not None

                if self.vb: print "Initializing",self.name,"model from scratch..."
                self.manual = False

            if self.flavor == 'Nebula':
                self.create_Nebula(position,SED)   # NB. Nebula position can be a list!
            else:
                self.create_Lens(position,SED)

        else:

            if self.vb: print "Initializing",self.name,"model from",template.name,"template..."

            if self.flavor == 'Nebula' and template.flavor == 'Nebula':
                assert self.K >= template.K
                self.spawn_Nebula(template)

            elif self.flavor == 'Lens' and template.flavor == 'Nebula':
                self.spawn_Lens(template)

        if self.vb:
            print "Initialization complete: "
            for component in self.srcs:
                print component
            print " "

        return None

# ----------------------------------------------------------------------------

    def create_Nebula(self,position,SED):

        # Start with a galaxy at the object centroid:
        if self.manual:
            galpos = position[0]
        else:
            galpos = position
        # Give it its fair share of the flux, and start faint:
        fudge = 0.2
        galSED = SED.copy() + 2.5*np.log10((self.K+1)/fudge)
        # Some standard shape and size parameters:
        re = resetre    # in arcsec, probably appropriate for the SQLS examples?
        q = 0.8     # axis ratio
        theta = 0.0 # degrees
        galshape = tractor.galaxy.GalaxyShape(re,q,theta)
        # Package up:
#        nebulousgalaxy = tractor.galaxy.ExpGalaxy(galpos,galSED,galshape)
#        nebulousgalaxy = tractor.galaxy.DevGalaxy(galpos,galSED,galshape)
        nebulousgalaxy = lenstractor.NebulousGalaxy(galpos,galSED,galshape)
        if self.vb: print nebulousgalaxy
        self.srcs.append(nebulousgalaxy)

        # Now add the point sources, with flux divided equally:
        for i in range(self.K):
            if self.manual:
                starpos = position[i+1]
            else:
                # Small random offsets from nebula centre:
                e = 0.2 # arcsec, MAGIC
                dx,dy = e*np.random.randn(2)/3600.0
                starpos = position.copy() + tractor.RaDecPos(dx,dy)
            starSED = SED.copy() + 2.5*np.log10((self.K+1)/fudge)
            # Package up:
            star = tractor.PointSource(starpos,starSED)
            if self.vb: print star
            self.srcs.append(star)

        return

# ----------------------------------------------------------------------------

    def spawn_Nebula(self,parent):

        self.srcs = []

        # Inherit the galaxy from the parent:
        nebulousgalaxy = parent.srcs[0]
        if self.vb: print nebulousgalaxy
        self.srcs.append(nebulousgalaxy)

        # Now inherit the point sources, and make more based on them, by
        # splitting the parent stars in 2 and separating them by a small
        # amount. Binary fission!
        parentstars = parent.srcs[1:]
        stars = []
        fluxratio = 0.2 # MAGIC
        k = 0
        while len(stars) < self.K:
            star1 = parentstars[k]
            parentBrightness = star1.getBrightness()
            star1.setBrightness(parentBrightness + 2.5*np.log10(1.0/(1.0-fluxratio)))
            stars.append(star1)
            if self.vb: print "Point source",star1
            k += 1
            if len(stars) < self.K:
                star2 = star1.copy()
                e = 0.1 # arcsec, MAGIC
                dx,dy = e*np.random.randn(2)/3600.0
                star2.setPosition(star2.getPosition() + tractor.RaDecPos(dx,dy))
                # BUG: not quite the right way to add small offsets in WCs...
                star2.setBrightness(parentBrightness + 2.5*np.log10(1.0/fluxratio))
                stars.append(star2)
                if self.vb: print "Point source",star2
        for star in stars:
            self.srcs.append(star)
        assert len(self.srcs) == (self.K + 1)

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
        re = resetre  # arcsec
        q = 0.8   # axis ratio
        theta = -phi # Note how mass/light misalignment is enabled.
        galshape = tractor.galaxy.GalaxyShape(re,q,theta)
        lensgalaxy = lenstractor.LensGalaxy(xd,md,galshape,thetaE,xshear)
        if self.vb: print lensgalaxy

        self.srcs.append(lenstractor.PointSourceLens(lensgalaxy, pointsource))

        return

# ----------------------------------------------------------------------------

    def spawn_Lens(self,parent):

        assert parent.flavor == 'Nebula'

        self.srcs = []

        # Inherit the lens galaxy from the parent Nebula:
        galaxy = parent.srcs[0]
        # Now inherit a point source, based on the Nebula's point sources!
        stars = parent.srcs[1:]

        # First, the lens galaxy light:
        xd = galaxy.getPosition()
        md = galaxy.getBrightness()
        galshape = galaxy.getShape()
        if self.vb:
            print "Galaxy's brightness = ",md[0]

        # NB this will fail if Nebula does not fit a galaxy well...


        # Now need thetaE (from image positions) and external shear.

        # First compute the image system centroid, unweighted:
        ra, dec = 0.0, 0.0
        for star in stars:
            radec = star.getPosition()
            ra += radec.ra
            dec += radec.dec
        ra, dec = ra/len(stars), dec/len(stars)
        centroid = tractor.RaDecPos(ra,dec)

        # Now we need the source centroid position *relative to the galaxy*:
        if self.vb:
            print "Solving for source position, shear etc..."
        self.deccorr = np.cos(xd.dec*deg2rad)
        # self.deccorr = 1.0
        self.xs0 = (ra - xd.ra)*self.deccorr*3600.0
        self.ys0 = (dec - xd.dec)*3600.0
        # First guess at source position is just the image centroid.
        xs = centroid
        if self.vb:
            print "Zeroth-order source displacement (arcsec) = ",self.xs0,", ",self.ys0

        # Make a rough initial estimate of Einstein radius:
        tE = 0.0
        for star in stars:
            tE += radec.distanceFrom(centroid)
        tE = tE*3600.0/len(stars)


        # - - - - - - - - - - - - - - - - - - - - - - - - -

        # Update source position according to Adri's math.
        # Remember to add back in xd ra and dec, and take care of dec!
        # (self.xs1, self.ys1) --> xs = tractor.RaDecPos(xd.ra + self.xs1/self.deccorr, xd.dec + self.ys1)

        # Set some parameters:

        ts = xd.distanceFrom(xs)*3600.0 # in arcsec
        magicdist=0.75*tE
        toodist=5*tE
        magicflat=0.3
        numb=len(stars)

        # PJM: Better would be to recurse until the system has converged.
        #      If it doesn't converge, there's a problem with the code that
        #      has to be fixed - so the code should *fail* in this
        #      eventuality.

        if self.vb:
            print "Number of images = ", numb
            print "Offset in Einstein radii = ",ts/tE
            print "Estimated Einstein Radius (arcsec) = ",tE

        # Do some basic checks to avoid obvious problems:

        if numb==2:
            # The system is a double:
            if ts >= toodist:
                print "Not a lens, the galaxy is way too far from the quasars!"
                assert False

            elif (ts >= magicdist or galshape.ab < magicflat or galshape.ab > 1./magicflat):
                if self.vb:
                    print "Galaxy may just be fitting noise, repositioning."

                # Reposition galaxy, reset axis ratio and p.a.:
                galpx, galpy, weight = 0.0, 0.0, 0.0
                pippo0 = stars[0].getBrightness()[0]
#                for star in stars:
                istar = 1
                while istar<len(stars):
                    star = stars[istar]
                    pippo1 = 10.0**(-0.4*star.getBrightness()[0]+0.4*pippo0)
                    galpx += star.getPosition().ra/pippo1
                    galpy += star.getPosition().dec/pippo1
                    weight += 1./pippo1
                    istar += 1
                galpx, galpy = galpx/weight, galpy/weight
                xd.ra, xd.dec = galpx, galpy

                # LT tries to reabsorb defects in the Nebula fit by playing
                # around with very elongated or localised galaxy components.
                # Reset these to ensure lens galaxy is roundish.
                galshape.ab = 1.1  # MAGIC reset axis ratio to be reasonably round...
                galshape.phi = 0.0 # MAGIC reset orientation to avoid mimicking the QSO residuals
                galshape.re = tE   # MAGIC reset eff.radius, in arcseconds, practically the extent of the img separation

                # Set zero shear, phi = p.a. of lens galaxy:
                gamma = 0.0
                gamma1, gamma2 = 0.0, 0.0
                phi = galshape.phi
                # Source position has already been set (xs)
                # Compute SIS total magnification for initial flux estimates
                mu = 2.0*tE/ts

            else:
                # Lens galaxy looks sensible, and has a well defined position -
                # so the shear cannot be zero. We have to solve for both the
                # shear and the Einstein radius, and in the process we get the
                # refined source position.

                xs,tE,gamma1,gamma2,mu = self.solve_for_initial_lens_mass_and_source_position(stars,xd,tE)


        elif numb==4:
            # The system is a quad - so we can solve for xd as well.

            if ts >= toodist:
                print "Not a lens, the galaxy is way too far from the quasars!"
                assert False

            elif (ts >= magicdist or galshape.ab < magicflat or galshape.ab > 1./magicflat):
                if self.vb:
                    print "Galaxy may just be fitting noise, repositioning."

                # Reposition galaxy using magnitudes in the first card, reset axis ratio and p.a.
                galpx, galpy, weight = 0.0, 0.0, 0.0
                pippo0 = stars[0].getBrightness()[0]
                for star in stars:
                    # pippo1 = 10.0**(-0.4*star.getBrightness()[0]+0.4*pippo0)
                    pippo1 = 1.
                    galpx += star.getPosition().ra/pippo1
                    galpy += star.getPosition().dec/pippo1
                    weight += 1./pippo1
                galpx, galpy = galpx/weight, galpy/weight
                xd.ra, xd.dec = galpx, galpy
                galshape.ab = 0.8 # magic reset to axis ratio
                galshape.phi = 0.0
                galshape.re = resetre

                # We need the centroid position *relative to the galaxy*, so we need to re-do this when resetting
                if self.vb:
                    print "Solving for source position, shear etc..."
                self.deccorr = np.cos(xd.dec*deg2rad)
                # self.deccorr = 1.0
                self.xs0 = (ra - xd.ra)*self.deccorr*3600.0
                self.ys0 = (dec - xd.dec)*3600.0
                self.xs0 *= -1.0 # To make x increase to the W, in a RH coord system, as in lens.py.
                # First guess at source position is just the image centroid.
                xs = centroid


            # Now, recursively optimize the centre of mass of the lens:

            ireccen=1
            Nreccen = 100 # No. of recursion steps for optimizing lens pos
            while ireccen<Nreccen:
                if self.vb:
                    print "Galaxy position = ", xd.ra,", ", xd.dec
                self.deccorr = np.cos(xd.dec*deg2rad)

                xs,tE,gamma1,gamma2,mu = self.solve_for_initial_lens_mass_and_source_position(stars,xd,tE)

                if self.vb:
                    print "tE (arcsec) = ", tE
                    print "total magnification = ", mu

                # Now we adjust the lens mass center:
                xd = self.solve_for_initial_lens_position(stars,xd,tE,gamma1,gamma2)

                ireccen += 1

        # Now we have everything we need to instantiate our Lens and Source models.


        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        # Instantiate ExternalShear object:

        gamma = (gamma1**2+gamma2**2)**0.5
        shreg = 0.0001
        phi = 0.5*np.arctan(gamma2/(gamma1+shreg)) # shreg helps in case gamma1==0, which is not expected to happen but one never knows..
        # phi = 0.5*np.arctan2(gamma2,(gamma1+shreg)) # shreg helps in case gamma1==0, which is not expected to happen but one never knows..
        # phi = 0.5*np.arctan2(gamma2,(gamma1+shreg)) # shreg helps in case gamma1==0, which is not expected to happen but one never knows..
        # PJM: tried arctan2. This code needs checking carefully against
        #      analytic results.
        phi = np.rad2deg(phi)
        xshear = lenstractor.ExternalShear(gamma,phi)
        # PJM: this is also potentially a source of error - a mismatch between your
        #      definition of phi and LensTractor's. Really, all this source estimation stuff
        #      should be a getLensedSource method in the LensGalaxy class (as I was writing).
        #      Not sure if that would help with the angle issues here, but it woudl keep
        #      the code modular...

        # Old workflow piece still valid
        ms = stars[0].getBrightness()
        # Start with just one point source's flux:
        pointsource = tractor.PointSource(xs,ms)
        # The Tractor likes starting with smaller fluxes...
        # mu = 10.0*mu
        # Add the other images' flux:
        for star in stars[1:]:
            pointsource.setBrightness(pointsource.getBrightness() + star.getBrightness())
        # Correct source brightness for approximate magnification:
        pointsource.setBrightness(pointsource.getBrightness() + 2.5*np.log10(mu))

        thetaE = lenstractor.EinsteinRadius(tE)
        if self.vb:
            print "Offset in Einstein radii = ",ts/tE
            print "Estimated Einstein Radius (arcsec) = ",thetaE
            print "Estimated shear amplitude gamma1 = ",gamma1
            print "Estimated shear amplitude gamma2 = ",gamma2
            print "Estimated shear amplitude gamma = ",gamma
            print "Estimated shear angle (degrees) = ",phi
        # Package into lensgalaxy:
        lensgalaxy = lenstractor.LensGalaxy(xd,md,galshape,thetaE,xshear)
        # if self.vb: print lensgalaxy
        # Note: this puts the lens mass where the galaxy light is!

        if self.vb: print pointsource

        self.srcs.append(lenstractor.PointSourceLens(lensgalaxy, pointsource))
        # assert False

        return

# ----------------------------------------------------------------------------

    def solve_for_initial_lens_mass_and_source_position(self,stars,xd,tE):

        #-- compute Ax etc
        Ax = Ay = sum0 = sum1 = sum2 = sum3 = sum4 = sum5 = A21 = A22 = A31 = A33 = B1 = B2 = 0.0
        #-- optimize for shear
        treg = 0.00001 # just not to have zero denominator in blabla/thetai
        for star in stars:
            # self.deccorr = 1.0
            stradec = star.getPosition()
            thetai = treg + stradec.distanceFrom(xd)*3600.0
            xi = (stradec.ra-xd.ra)*self.deccorr*3600.0
            yi = (stradec.dec-xd.dec)*3600.0
            xi *= -1.0 # To make x increase to the W, in a RH coord system, as in lens.py.
            Ax += xi/thetai # Units should be the same above and below
            Ay += yi/thetai # And here too
            sum0 += thetai
            sum1 += (xi**2 - yi**2)/thetai
            sum2 += thetai**2
            sum3 += yi**2 - xi**2
            sum4 += xi*yi/thetai
            sum5 += xi*yi
        #-- params of tE as a function of gamma
        A11 = len(stars) -(Ax**2)/len(stars) -(Ay**2)/len(stars)
        A12 = sum1 - Ax*self.xs0 +Ay*self.ys0
        A13 = 2.0*sum4 -Ax*self.ys0 - Ay*self.xs0
        B1 = sum0 - Ax*self.xs0 - Ay*self.ys0
        A21 = (self.ys0*Ay-self.xs0*Ax+sum1)/len(stars) # self.xs0, self.ys0 are in arcseconds
        A22 = sum2/len(stars) - self.xs0**2 - self.ys0**2
        B2 = self.ys0**2 - self.xs0**2 -sum3/len(stars)
        gamma1 = -(B2-A21*tE)/A22
        A31 = -(Ax*self.ys0+Ay*self.xs0-2*sum4)/len(stars)
        A33 = A22
        B3 = +2.0*sum5/len(stars) -2.0*self.xs0*self.ys0/len(stars)
        gamma2 = -(B3 - A31*tE)/A33
        #-- solve recursively in tE, gamma1, gamma2
        irec = 1
        NrectE = 200  # No. of recursion steps for optimizing thetaE
        while irec<NrectE:
            tE = (B1 + gamma1*A12 + A13*gamma2)/A11
            gamma1 = -(B2 - A21*tE)/A22
            gamma2 = -(B3 - A31*tE)/A33
            irec += 1

        # Now compute the source position, and magnification:

        self.xs1 = self.ys1 = 0.0
        mu = 0.0 # for the SIS+XS total magnification
        muinvreg = 0.001 # to cope with possibly infinite magnification (although vey unlikely)
        for star in stars:
            # self.deccorr = 1.0
            stradec = star.getPosition()
            thetai = treg + stradec.distanceFrom(xd)*3600.0
            dxs = (stradec.ra-xd.ra)*self.deccorr*3600.0
            dys = (stradec.dec-xd.dec)*3600.0
            dxs *= -1.0 # To make x increase to the W, in a RH coord system, as in lens.py.
            self.xs1 += dxs*(1.0 - tE/thetai + gamma1) + gamma2*dys # lens equation
            self.ys1 += dys*(1.0 - tE/thetai - gamma1) + gamma2*dxs # lens equation
            muinv = 1.0 - gamma1**2 - gamma2**2 +(tE/thetai)*(-1.0 - gamma1*(dxs**2+dys**2)/thetai - 2.0*gamma2*dxs*dys/thetai**2)
            mu += 1.0/(np.abs(muinv) + muinvreg)
        self.xs1, self.ys1 = self.xs1/len(stars), self.ys1/len(stars) # Nb. in arcsec
        ras1, decs1 = -self.xs1/(self.deccorr*3600.0), self.ys1/3600.0 # Nb. in degrees, on sky

        xs = tractor.RaDecPos(xd.ra + ras1, xd.dec + decs1)

        return xs,tE,gamma1,gamma2,mu

# ----------------------------------------------------------------------------

    def solve_for_initial_lens_position(self,stars,xd,tE,gamma1,gamma2):

        learnrate = 0.1
        d11 = d12 = 0.0
        d21 = d22 = 0.0
        dchi2xd = dchi2yd = 0.0
        squares = 0.0 # used for source-plane chi2, i.e. variance in pred.source positions
        treg = 0.00001 # just not to have zero denominator in blabla/thetai
        for star in stars:
            stradec = star.getPosition()
            dxs = (stradec.ra-xd.ra)*self.deccorr*3600.0
            dys = (stradec.dec-xd.dec)*3600.0
            dxs *= -1.0 # To make x increase to the W, in a RH coord system, as in lens.py.
            thetai = treg + (dxs**2 + dys**2)**0.5
            # thetai = stradec.distanceFrom(xd) # in degrees
            # thetai = treg+ thetai*3600.0 # in arcseconds
            d11 = -1.0 - gamma1 + tE*(dys**2)/(thetai**3)
            d12 = -gamma2 - tE*dys*dxs/(thetai**3)
            d21 = d12
            d22 = -1.0 + gamma1 + tE*(dxs**2)/(thetai**3)
            xsi = dxs*(1.0 -tE/thetai + gamma1) + gamma2*dys # in arcseconds. Note gamma sign convention...
            ysi = dys*(1.0 -tE/thetai - gamma1) + gamma2*dxs # in arcseconds. Note gamma sign convention...
            dchi2xd += (2./len(stars))*(xsi*d11 +ysi*d12 - self.xs1*d11 - self.ys1*d12) # in arcseconds, self.xs1 ans self.ys1 are in degrees
            dchi2yd += (2./len(stars))*(xsi*d21 +ysi*d22 - self.xs1*d21 - self.ys1*d22) # in arcseconds
            # squares += (1./len(stars))*(xsi**2 + ysi**2)
        dispx = dchi2xd/self.deccorr # in arcseconds
        dispy = dchi2yd # in arcseconds
        dispx *= -1.0 # To make RA increase to the E, in a LH coord system, as in lens.py.
        displamp = (dispx**2 + dispy**2)**0.5 # in arcseconds
        displamp = displamp/tE # pure number
        if self.vb:
            # print "source-plane chi2  = ",squares - (self.xs1**2 +self.ys1**2)
            print "tE (arcsec) = ", tE
            print "displacement/tE = ", displamp
            print "deflector's adjustments (units of tE) = ",dispx*learnrate/((learnrate+displamp)*tE),dispy*learnrate/((learnrate+displamp)*tE)
        xd.ra, xd.dec = xd.ra -(dispx/3600.0)*learnrate/(learnrate+displamp), xd.dec -(dispy/3600.0)*learnrate/(learnrate+displamp) # in degrees

        return xd

# ----------------------------------------------------------------------------

    def get_positions_from(self,catalog,SED):

        positions = []
        Nbands = SED.numberOfParams()

        # Open up LT output catalog format file and read positions, assuming
        # hard-coded column numbers from the PS1 example:
        # Need to change that for other surveys!
        x = np.loadtxt(catalog)

        # Init catalogs are always for Nebulae models.
        # Parameters are something like:
        #
        #  1        catalog.source0.pos.ra   1.478442646221303107e+02
        #  2       catalog.source0.pos.dec   2.658710567091902277e+01
        #  3  catalog.source0.brightness.g   1.838935729060930768e+01
        #  4  catalog.source0.brightness.i   1.781994514057002377e+01
        #  5  catalog.source0.brightness.r   1.791359861295001821e+01
        #  6  catalog.source0.brightness.u   1.861519625801571465e+01
        #  7  catalog.source0.brightness.z   1.790694668268868384e+01
        #  8      catalog.source0.shape.re   2.356612496948709001e-01
        #  9      catalog.source0.shape.ab   2.286815004837887244e+00
        # 10     catalog.source0.shape.phi  -4.142924582249651877e+01
        # 11        catalog.source1.pos.ra   1.478440593629958073e+02
        # 12       catalog.source1.pos.dec   2.658722018688290234e+01
        # 13  catalog.source1.brightness.g   1.889259658866843949e+01
        # 14  catalog.source1.brightness.i   1.835678754819763014e+01
        # 15  catalog.source1.brightness.r   1.844249717329058669e+01
        # 16  catalog.source1.brightness.u   1.896585704168741060e+01
        # 17  catalog.source1.brightness.z   1.838253660524710398e+01
        # 18        catalog.source2.pos.ra   1.478440533891872235e+02
        # 19       catalog.source2.pos.dec   2.658722428423748241e+01
        # 20  catalog.source2.brightness.g   1.885674019434759785e+01
        # 21  catalog.source2.brightness.i   1.835242649962165373e+01
        # 22  catalog.source2.brightness.r   1.844604212419161726e+01
        # 23  catalog.source2.brightness.u   1.893952465584669653e+01
        # 24  catalog.source2.brightness.z   1.829995958971434789e+01

        # ie Npars = (5 + Nbands) + K*(2 + Nbands)
        # where K is the no. of point sources in the Nebula.

        numerator = (len(x) - (5 + Nbands))
        denominator = (2 + Nbands)
        assert numerator%denominator == 0
        Npos = numerator/denominator
        if (Npos < 1 or Npos > 4):
            print "ERROR: unrecognised Nebula model in catalog. K = ",Npos
            assert False

        # Set up Nebula model:
        self.K = Npos

        # Create position objects:
        j = 0
        positions.append(tractor.RaDecPos(x[j],x[j+1]))
        j += 3
        for i in range(Npos):
            j += 2 + Nbands
            positions.append(tractor.RaDecPos(x[j],x[j+1]))
        assert len(positions) == (self.K+1)

        return positions

# ----------------------------------------------------------------------------

    def plot(self,wcs,band):

        if self.flavor == 'Nebula':
            self.plot_Nebula(wcs,band)
        else:
            self.plot_Lens(wcs,band)

        return

# ----------------------------------------------------------------------------

    def plot_Nebula(self,wcs,band):

        galaxy = self.srcs[0]
        stars = self.srcs[1:]

        # Plot galaxy as orange circle:
        radec = galaxy.getPosition()
        x,y = wcs.positionToPixel(radec)
        SED = galaxy.getBrightness()
        plotmag = (20.0 - SED.getMag(band))*50 # MAGIC 20,50
        plt.scatter(x,y,color='orange',s=plotmag,alpha=0.3)

        # Plot point sources as cyan circles:
        for star in stars:
            radec = star.getPosition()
            x,y = wcs.positionToPixel(radec)
            SED = star.getBrightness()
            plotmag = (20.0 - SED.getMag(band))*50 # MAGIC 20,50
            plt.scatter(x,y,color='cyan',s=plotmag,alpha=0.3)

        return

# ----------------------------------------------------------------------------

    def plot_Lens(self,wcs,band):

        # Plot lens galaxy as orange circle:
        # lensgalaxy = self.srcs[0]
        # source = self.srcs[1:]

        # Get image positions:
        # FUNCTIONCALL

        # Plot point images as cyan circles:
        # for ptimg in ptimgs:



        return

# ============================================================================

if __name__ == '__main__':

   pass
