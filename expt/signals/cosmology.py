from numpy import log10
import _cosmo

__all__ = ['Cosmology', 'StdCdlCosmology']

# TODO:  Implement this as an extension object to avoid the set_frw
# calls.

class Cosmology(object):

    def __init__(self, h, O_m, O_lambda):
        self.h = h
        self.O_m = O_m
        self.O_lambda = O_lambda
        self.O_k = 1. - O_m - O_lambda

    def ldist(self, z):
        """
        Return the luminosity distance to redshift z, in Mpc.
        """
        _cosmo.set_frw_cosmo(self.h, self.O_m, self.O_lambda)
        try:
            len(z)
            ldvals = _cosmo.nldist(z)
            return ldvals / 1.e5
        except TypeError:
            return _cosmo.ldist(z) / 1.e5

    def mu(self, z):
        """
        Return the distance modulus to redshift z.
        This will work for a vector of z values, as well as a scalar.
        """
        _cosmo.set_frw_cosmo(self.h, self.O_m, self.O_lambda)
        try:
            len(z)
            ldvals = _cosmo.nldist(z)
            return 5*log10(ldvals)
        except TypeError:
            return _cosmo.mu_z(z)

    def vol_elem(self, z):
        """
        Return the volume per unit redshift per unit steradian.
        """
        _cosmo.set_frw_cosmo(self.h, self.O_m, self.O_lambda)
        return _cosmo.vol_elem(z)

    def dlbt(self, z):
        """
        Return the dimensionless look-back time to redshift z
        (in units of 1/H_0).
        """
        _cosmo.set_frw_cosmo(self.h, self.O_m, self.O_lambda)
        return _cosmo.dlbt_z(z)

    def dage(self):
        """
        Return the dimensionless age of the universe for the 
        current cosmology (in units of 1/H_0).
        """
        _cosmo.set_frw_cosmo(self.h, self.O_m, self.O_lambda)
        return _cosmo.dage()


class StdCdlCosmology(Cosmology):

    def __init__(self, h, O_m, O_lambda, F_fid):
        self.F_fid = F_fid
        _cosmo.set_ffid(F_fid)
        Cosmology.__init__(self, h, O_m, O_lambda)

    def set_lum(self, lum):
        """
        Set the dimensionless standard-candle luminosity.
        """
        self.lum = lum
        _cosmo.set_lum(lum)

    def flux2z(self, F):
        """
        Return the redshift such that a standard-candle source
        with have flux F.
        """
        _cosmo.set_ffid(self.F_fid)
        _cosmo.set_lum(self.lum)
        return _cosmo.flux2z(F)