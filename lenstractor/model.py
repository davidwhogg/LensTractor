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
                
        # Start with an De Vaucouleurs galaxy at the object centroid:
        if self.manual:
            galpos = position[0]
        else:
            galpos = position
        # Give it its fair share of the flux, and start faint:
        fudge = 0.2 # MAGIC
        galSED = SED.copy() + 2.5*np.log10((self.K+1)/fudge)
        # Some standard shape and size parameters:
        re = 0.5    # arcsec
        q = 0.8     # axis ratio
        theta = 0.0 # degrees
        galshape = tractor.sdss_galaxy.GalaxyShape(re,q,theta)
        # Package up:
        nebulousgalaxy = tractor.sdss_galaxy.DevGalaxy(galpos,galSED,galshape)
        if self.vb: print nebulousgalaxy
        self.srcs.append(nebulousgalaxy)

        # Now add the point sources, with flux divided equally:
        for i in range(self.K):
            if self.manual:
                starpos = position[i+1]
            else:
                # Small random offsets from nebula centre:
                e = 0.1 # arcsec, MAGIC
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
        
        # Inherit the De Vaucouleurs galaxy from the parent:
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
        re = 0.5  # arcsec
        q = 0.8   # axis ratio
        theta = -phi # Note how mass/light misalignment is enabled.
        galshape = tractor.sdss_galaxy.GalaxyShape(re,q,theta)
        lensgalaxy = lenstractor.LensGalaxy(xd,md,galshape,thetaE,xshear)
        if self.vb: print lensgalaxy

        self.srcs.append(lenstractor.PointSourceLens(lensgalaxy, pointsource))

        return
        
# ----------------------------------------------------------------------------
    
    def spawn_Lens(self,parent):
        
        assert parent.flavor == 'Nebula'
        
        self.srcs = []
        
        # Inherit the De Vaucouleurs lens galaxy from the parent Nebula:
        galaxy = parent.srcs[0]
        # Now inherit a point source, based on the Nebula's point sources!
        stars = parent.srcs[1:]
        
        # First, the lens galaxy light:
        xd = galaxy.getPosition()
        md = galaxy.getBrightness()
        galshape = galaxy.getShape()
        if self.vb:
            print "Galaxy's brightness = ",md[0]
#        
        # Now need thetaE (from image positions) and external shear. 
        # Compute the image system centroid and, from that, the Einstein radius:
        ra, dec = 0.0, 0.0
        for star in stars: 
            radec = star.getPosition()
            ra += radec.ra
            dec += radec.dec
        ra, dec = ra/len(stars), dec/len(stars)
        centroid = tractor.RaDecPos(ra,dec)
#       for the new stuff: we need centroid position relative to galaxy
        xs0, ys0 = ra-xd.ra, dec-xd.dec # there are in degrees!
        if self.vb:
            print "Zeroth-order source displacement = ",xs0,", ",ys0
#       Old spawn
        tE = 0.0
        for star in stars: 
            tE += radec.distanceFrom(centroid)
        tE = tE*3600.0/len(stars)
#       Old spawn        
#        # Associate the lens light ellipticity with the lens potential shear...
#        # (Warning, this could be crazy!)
#        q = galshape.ab
#        gamma = 0.2#*(1-q)/(1+q) # MAGIC 0.5
#        phi   = galshape.phi # deg
#        xshear = lenstractor.ExternalShear(gamma,phi)
#        tE = tE*(1-gamma)*0.5 # start with separation sensibly smaller than the naively estimated one
        
        # First guess at source position is just the number-averaged centroid        
        # will need to update it (xs1, ys1)--> xs = tractor.RaDecPos(xs1+xd.ra,ys1+xd.dec)                 
        xs = centroid
        ts = xd.distanceFrom(centroid)*3600.0
#       New workflow
#        irep=1
#        while irep<10:
#        #indent!
        magicdist=0.5*tE
        numb=len(stars)
        if self.vb:
            print "Number of peaks = ", numb
        if numb==2: # spawn double
            if ts >= magicdist:
                #-- reposition galaxy, reset axis ratio and p.a.
                galpx, galpy, weight = 0.0, 0.0, 0.0
                for star in stars:
                    galpx += star.getPosition().ra/star.getBrightness()[0]
                    galpy += star.getPosition().dec/star.getBrightness()[0]
                    weight += 1./star.getBrightness()[0]
                galpx, galpy = galpx/weight, galpy/weight
                xd.ra, xd.dec = galpx, galpy
                galshape.ab = 0.8 # magic reset to axis ratio
                galshape.phi = 0.0                    
                #-- set zero shear, phi = p.a.
                gamma = 0.0
                phi = galshape.phi
                #-- leave xs as it is
            else:
                #-- compute Ax etc
                Ax = Ay = sum1 = sum2 = sum3 = sum4 = sum5 = A21 = A22 = A31 = A33 = B1 = B2 = 0.0
                #-- optimize for shear
                treg = 0.001 # just not to have zero denominator in blabla/thetai
                for star in stars:
                    stradec = star.getPosition()
                    thetai = treg + stradec.distanceFrom(xd)*3600
                    xi = (stradec.ra-xd.ra)*3600
                    yi = (stradec.dec-xd.dec)*3600
                    Ax += xi/thetai # units should be the same above and below
                    Ay += yi/thetai # same here
                    sum1 += (xi**2 - yi**2)/thetai
                    sum2 += thetai**2
                    sum3 += yi**2 - xi**2
                    sum4 += xi*yi/thetai
                    sum5 += xi*yi
#                Ax, Ay, sum1, sum2, sum3, sum4, sum5 = Ax/len(stars), Ay/len(stars), sum1/len(stars), sum2/len(stars), sum3/len(stars), sum4/len(stars), sum5/len(stars)
                A21 = (ys0*3600*Ay-xs0*3600*Ax+sum1)/len(stars) # I'm converting xs0, ys0 in arcseconds like the other quantities
                A22 = sum2/len(stars) -(xs0*3600)**2 -(ys0*3600)**2
                B2 = (ys0*3600)**2 - (xs0*3600)**2 +sum3/len(stars)
                gamma1 = (B2-A21*tE)/A22
                A31 = (Ax*ys0*3600-Ay*xs0*3600-2*sum4)/len(stars)
                A33 = 2*sum2/len(stars) - (ys0*3600)**2 - (xs0*3600)**2
                B3 = 2*sum5/len(stars) -xs0*ys0*(3600**2)
                gamma2 = (B3 - A31*tE)/A33
                #-- compute gamma, phi, xs1, ys1, update xs
                shreg = 0.0001
                gamma = (gamma1**2+gamma2**2)**0.5
                phi = 0.5*np.arctan(gamma2/(gamma1+shreg)) # shreg helps in case gamma1==0, which is not expected to happen but one never knows..
                xs1 = ys1 = 0.0
                for star in stars:
                    stradec = star.getPosition()
                    thetai = treg + stradec.distanceFrom(xd)*3600
                    xs1 += (stradec.ra-xd.ra)*(1.0-tE/thetai-gamma1) - gamma2*(stradec.dec-xd.dec) # lens equation
                    ys1 += (stradec.dec-xd.dec)*(1.0-tE/thetai+gamma1) - gamma2*(stradec.ra-xd.ra) # lens equation
                xs1, ys1 = xs1/len(stars), ys1/len(stars)
                xs = tractor.RaDecPos(xs1+xd.ra,ys1+xd.dec) 
        else: # spawn quad, practically same as above
            if ts >= magicdist:
                if self.vb:
                    print "Galaxy may be fictitious, repositioning."
                #-- reposition galaxy using magnitudes in the first card, reset axis ratio and p.a.
                galpx, galpy, weight = 0.0, 0.0, 0.0
                pippo0 = stars[0].getBrightness()[0]
                for star in stars:
                    pippo1 = 10**(-0.4*star.getBrightness()[0]+0.4*pippo0)
                    galpx += star.getPosition().ra/pippo1
                    galpy += star.getPosition().dec/pippo1
                    weight += 1./pippo1
                galpx, galpy = galpx/weight, galpy/weight
                xd.ra, xd.dec = galpx, galpy
                galshape.ab = 0.8 # magic reset to axis ratio
                galshape.phi = 0.0                    
            #-- else: do nothing
        if self.vb:
            print "Galaxy position = ", xd.ra,", ", xd.dec
            #-- now compute Ax etc
            Ax = Ay = sum1 = sum2 = sum3 = sum4 = sum5 = A21 = A22 = A31 = A33 = B1 = B2 = 0.0
            #-- optimize for shear
            treg = 0.001 # just not to have zero denominator in blabla/thetai
            for star in stars:
                stradec = star.getPosition()
                thetai = treg + stradec.distanceFrom(xd)*3600
                xi = (stradec.ra-xd.ra)*3600
                yi = (stradec.dec-xd.dec)*3600
                Ax += xi/thetai # units should be the same above and below
                Ay += yi/thetai # same here
                sum1 += (xi**2 - yi**2)/thetai
                sum2 += thetai**2
                sum3 += yi**2 - xi**2
                sum4 += xi*yi/thetai
                sum5 += xi*yi
#            Ax, Ay, sum1, sum2, sum3, sum4, sum5 = Ax/len(stars), Ay/len(stars), sum1/len(stars), sum2/len(stars), sum3/len(stars), sum4/len(stars), sum5/len(stars)
#            if self.vb:
#                print "Shear optimization parameters = ", Ax,", ", Ay,", ", sum1,", ", sum2,", ", sum3,", ", sum4,", ", sum5
#
            A21 = (ys0*3600*Ay-xs0*3600*Ax+sum1)/len(stars) # I'm converting xs0, ys0 in arcseconds like the other quantities
            A22 = sum2/len(stars) -(xs0*3600)**2 -(ys0*3600)**2
            B2 = (ys0*3600)**2 - (xs0*3600)**2 +sum3/len(stars)
            gamma1 = (B2-A21*tE)/A22
            A31 = (Ax*ys0*3600-Ay*xs0*3600-2*sum4)/len(stars)
            A33 = 2*sum2/len(stars) - (ys0*3600)**2 - (xs0*3600)**2
            B3 = 2*sum5/len(stars) -xs0*ys0*(3600**2)
            gamma2 = (B3 - A31*tE)/A33
            #-- compute gamma, phi, update xs
            shreg = 0.0001
            gamma = (gamma1**2+gamma2**2)**0.5
            phi = 0.5*np.arctan(gamma2/(gamma1+shreg)) # shreg helps in case gamma1==0, which is not expected to happen but one never knows..
            xs1 = ys1 = 0.0
            for star in stars:
                stradec = star.getPosition()
                thetai = treg + stradec.distanceFrom(xd)*3600
                xs1 += (stradec.ra-xd.ra)*(1.0-tE/thetai-gamma1) - gamma2*(stradec.dec-xd.dec) # lens equation
                ys1 += (stradec.dec-xd.dec)*(1.0-tE/thetai+gamma1) - gamma2*(stradec.ra-xd.ra) # lens equation
            xs1, ys1 = xs1/len(stars), ys1/len(stars)
            xs = tractor.RaDecPos(xs1+xd.ra,ys1+xd.dec)
#
#           irep++
# Iterate? I.e. optimize in tE given gamma1, gamma2, re-solve in gamma1, gamma2 etc.
#
# instantiate ExternalShear object
        phi = np.rad2deg(phi)
        xshear = lenstractor.ExternalShear(gamma,phi)
        
#       Old workflow piece still valid
        ms = stars[0].getBrightness()
        # Start with just one point source's flux:
        pointsource = tractor.PointSource(xs,ms)
        # Add the other images' flux:
        for star in stars[1:]: 
            pointsource.setBrightness(pointsource.getBrightness() + star.getBrightness())
        # Correct source brightness for approximate magnification:
        mu = 2*tE/ts
        pointsource.setBrightness(pointsource.getBrightness() + 2.5*np.log10(mu))
#
        thetaE = lenstractor.EinsteinRadius(tE)
        if self.vb:
            print "Offset in Enstein radii = ",ts/tE
            print "Estimated Einstein Radius (arcsec) = ",thetaE
            print "Estimated shear amplitude gamma1 = ",gamma1
            print "Estimated shear amplitude gamma2 = ",gamma2
            print "Estimated shear amplitude gamma = ",gamma
            print "Estimated shear angle (degrees) = ",phi
        # Package into lensgalaxy:
        lensgalaxy = lenstractor.LensGalaxy(xd,md,galshape,thetaE,xshear)
#        if self.vb: print lensgalaxy
        # Note: this puts the lens mass where the galaxy light is!
        

        if self.vb: print pointsource

        self.srcs.append(lenstractor.PointSourceLens(lensgalaxy, pointsource))
        return
        
# ----------------------------------------------------------------------------
    
    def get_positions_from(self,catalog,SED):
    
        positions = []
        Nbands = SED.numberOfParams()
        
        # Open up LT output catalog format file and read positions, assuming 
        # hard-coded column numbers:
        x = np.loadtxt(catalog)
        
        if len(x) == 23:
            Npos = 4
            self.K = Npos
        elif len(x) == 13:
            Npos = 2
            self.K = Npos
        else:
            print "ERROR: unrecognised Nebula model in catalog, with ",len(x)," parameters"
            assert False
        
        # Catalog must match this model!
        assert Npos == self.K #overridden by self.K initialisation above
        
        # Create position objects:
        j = 0
        positions.append(tractor.RaDecPos(x[j],x[j+1]))
        j += 3
        for i in range(Npos):
            j += 2 + Nbands
            positions.append(tractor.RaDecPos(x[j],x[j+1]))

        return positions
        
# ============================================================================

if __name__ == '__main__':

   pass
