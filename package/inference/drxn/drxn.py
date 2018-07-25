from numpy import array, zeros, dot, sin, cos, log, pi, sqrt
import fisher_coinc as fc

__all__ = ['Direction', 'FDirection',
           'doublet_lbf', 'triplet_lbf', 'multiplet_lbf',
           'fisher_lml', 'fisher_fit']

class Direction(object):
    """A direction to a point on the unit sphere."""

    def __init__(self, lng=None, lat=None, ra=None, dec=None, deg=False,
                 uvec=None):
        """
        Initialize a direction given its angular coordinates.

        lng, lat = longitude and latitude
        ra, dec = right acension and declination
        deg = Units flag; False->radians, True->degrees

        uvec = unit vector specifying direction

        Only one of (long,lat) or (ra,dec) or uvec may be provided.
        """
        # Store the direction and components of a unit vector in
        # the coordinate system of the input angles.
        if uvec is not None:
            if lng is not None or lat is not None or ra is not None \
               or dec is not None:  # Don't test deg so subclasses can use it
                raise ValueError, 'Specify only one type of direction!'
            self.uvec = uvec
            self.long, self.lat = None, None
            self.ra, self.dec = None, None
        elif lng is not None and lat is not None:
            if ra is not None or dec is not None:
                raise ValueError, 'Use only long/lat OR ra/dec!'
            if deg:
                lng, lat = pi*lng/180., pi*lat/180.
            self.long, self.lat = lng, lat
            self.ra, self.dec = None, None
            self.uvec = zeros(3,float)
            clat = cos(lat)
            self.uvec[0] = cos(lng) * clat
            self.uvec[1] = sin(lng) * clat
            self.uvec[2] = sin(lat)
        elif ra is not None and dec is not None:
            if lng is not None or lat is not None:
                raise ValueError, 'Use only long/lat OR ra/dec!'
            if deg:
                ra, dec = pi*ra/180., pi*dec/180.
            self.ra, self.dec = ra, dec
            self.long, self.lat = None, None
            self.uvec = zeros(3,float)
            cdec = cos(dec)
            self.uvec[0] = cos(ra) * cdec
            self.uvec[1] = sin(ra) * cdec
            self.uvec[2] = sin(dec)
        else:
            print lng, lat, ra, dec
            raise ValueError, 'Illegal set of arguments!'

    def __mul__(self, drxn):
        """The dot product between direction unit vectors."""
        return sum(self.uvec*drxn.uvec)

    def __str__(self):
        s = ''
        if self.long is not None:
            s += 'long, lat (rad): %g %g\n' % (self.long, self.lat)
        elif self.ra is not None:
            s = 'RA, dec (rad): %g %g\n' % (self.ra, self.dec)
        return s + str(self.uvec)


class FDirection(Direction):
    """A "Fisher direction" object, storing the direction to a point on
    the unit sphere, and a description of the direction uncertainty
    in terms of the Fisher distribution, a spherical generalization
    of the normal distribution."""

    log4pi = log(4*pi)

    def __init__(self, sig=None, lng=None, lat=None, ra=None, dec=None, deg=False,
                 uvec=None, kappa=None):
        """
        Initialize a Fisher direction given its angular uncertainty
        and direction.

        sig = polar anglular radius of 68.3% uncertainty region
        kappa = Fisher distribution concentration parameter
        Only one of sig or kappa may be provided.
        
        lng, lat = longitude and latitude
        ra, dec = right acension and declination
        deg = Units flag; False->radians, True->degrees

        uvec = Unit vector specifying direction

        Only one of (long,lat) or (ra,dec) or uvec may be provided.
        """
        Direction.__init__(self, lng, lat, ra, dec, deg, uvec)
        if kappa:
            self.kappa = kappa
            self.sig = None
        else:
            if deg: sig = pi*sig/180.
            self.sig = sig
            self.kappa = fc.sig_to_kappa(sig)
        self.lsk = fc.lsinhc(self.kappa)
        self.const = self.lsk + self.log4pi

    def llike(self, drxn=None, lng=None, lat=None, ra=None, dec=None, uvec=None):
        """The log-likelihood that the given direction is the true direction
        for this FDrxn instance."""
        if drxn:
            if lng or lat or ra or dec or uvec:
                raise ValueError, 'Use only drxn or long/lat OR ra/dec OR uvec!'
            uvec = drxn.uvec
        if lng is not None and lat is not None:
            if ra is not None or dec is not None or uvec is not None:
                raise ValueError, 'Use only long/lat OR ra/dec OR uvec!'
            uvec = zeros(3,float)
            clat = cos(lat)
            uvec[0] = cos(lng) * clat
            uvec[1] = sin(lng) * clat
            uvec[2] = sin(lat)
        elif ra is not None and dec is not None:
            if lng is not None or lat is not None or uvec is not None:
                raise ValueError, 'Use only long/lat OR ra/dec OR uvec!'
            uvec = zeros(3,float)
            cdec = cos(dec)
            uvec[0] = cos(ra) * cdec
            uvec[1] = sin(ra) * cdec
            uvec[2] = sin(dec)
        elif uvec:
            if long or lat or ra or dec:
                raise ValueError, 'Use only long/lat OR ra/dec OR uvec!'
        return self.kappa*dot(uvec,self.uvec) - self.const

    def __str__(self):
        s = Direction.__str__(self) + '\n'
        s += 'kappa = %g' % self.kappa
        return s

    def __mul__(self, drxn):
        """The dot product between direction unit vectors."""
        return sum(self.uvec*drxn.uvec)


def doublet_lbf(drxn1, drxn2):
    """Log Bayes factor favoring coincidence of 2 directions."""
    return fc.doublet_lbf(drxn1.uvec, drxn1.kappa, drxn1.lsk,
                          drxn2.uvec, drxn2.kappa, drxn2.lsk)

def triplet_lbf(drxn1, drxn2, drxn3):
    """Log Bayes factor favoring coincidence of 3 directions."""
    return fc.triplet_lbf(drxn1.uvec, drxn1.kappa, drxn1.lsk,
                          drxn2.uvec, drxn2.kappa, drxn2.lsk,
                          drxn3.uvec, drxn3.kappa, drxn3.lsk)

def multiplet_lbf(drxns):
    """Log Bayes factor favoring coincidence of a sequence
    of directions."""
    nd = len(drxns)
    uvecs = array([drxn.uvec for drxn in drxns])
    uvecs = uvecs.transpose() # For Fortran compatibility
    kappas = array([drxn.kappa for drxn in drxns])
    lsks = array([drxn.lsk for drxn in drxns])
    return fc.multiplet_lbf(uvecs, kappas, lsks)

def fisher_lml(drxns):
    """Log maximum likelihood ratio for fit of Fisher directions to a common
    direction.  This is equivalent to the NVO XMatch algorithm."""
    # Collect kappa*uvec
    kuvecs = array([drxn.kappa*drxn.uvec for drxn in drxns])
    # Add all the scaled unit vectors
    kuvec = sum(kuvecs,0)
    ksum = sum(array([drxn.kappa for drxn in drxns]))
    return ksum - sqrt(sum(kuvec**2))
    
def fisher_fit(drxns):
    """
    Maximum likelihood fit of a sequence of Fisher directions to a common
    direction.  This is equivalent to the NVO XMatch algorithm.

    lml, drxn = fisher_fit(drxns)

    lmll = log maximum likelihood ratio for the fit
    drxn = best-fit Fisher direction.
    """
    # Collect kappa*uvec
    kuvecs = array([drxn.kappa*drxn.uvec for drxn in drxns])
    # Add all the scaled unit vectors
    kuvec = sum(kuvecs,0)
    kappa = sqrt(sum(kuvec**2))
    ksum = sum(array([drxn.kappa for drxn in drxns]))
    return ksum - kappa, FDirection(kappa=kappa, uvec=kuvec/kappa)
    
