from math import exp
from roots import zbrac, zbrent

class pdgrid1:
    """1D probability density grid."""

    def __init__(self, lo, hi, n, step="linear"):
        """Initialize a grid from lo to hi with n steps,
        either linearly or logarithmically spaced."""
        pass


def qgt1(a, crit):
    """Quadrature Greater Than (1d)
    Return the trapezoid rule quadrature of the contents of 1d array a
    for the region where the integrand is above the value crit.
    The array and crit contain the *log* of integrand values.
    Step size factors are *not* included and should be provided by
    the calling routine."""

    if len(a.shape) != 1:
        raise ValueError, '1D array required!'
    q = 0.
    amaxp = a[0]    # Hereafter, the previous max
    amax = amaxp
    for i in range(a.shape[0]-1):
        if a[i] > crit and a[i+1] > crit:
            amax = max(amaxp, a[i], a[i+1])
            q = q * exp(amaxp-amax) + \
                0.5 * (exp(a[i]-amax) + exp(a[i+1]-amax))
            amaxp = amax
            # print 'both: ',i, a[i], a[i+1], amax, q
        elif a[i] > crit or a[i+1] > crit:
            hi = max(a[i], a[i+1])
            lo = min(a[i], a[i+1])
            amax = max(amaxp, hi)
            lo, hi, c = exp(lo-amax), exp(hi-amax), exp(crit-amax)
            q = q * exp(amaxp-amax) + \
                0.5 * (hi + c) * (hi - c) / (hi - lo)
            # print 'one: ',i, a[i], a[i+1], amax, q
    return amax, q

from numpy import *
d = 0.01
xvals = arange(-5.,5.,d)
a = []
for x in xvals:
    a.append( -0.5*x**2 )
a = array(a) - log(sqrt(2.*pi))

crit = -2. - log(sqrt(2.*pi))
am, q = qgt1(a, crit)
print am, q, d*q*exp(am)
