from pylab import *
from numpy import *
from cosmology import *

zvals = linspace(0., 10., 6)
print 'z vals:', zvals
print

cosmo = Cosmology(0.7, 0.3, 0.7)
print cosmo.ldist(2), cosmo.ldist(6)
print cosmo.ldist(zvals)
print cosmo.vol_elem(1.), cosmo.dlbt(1.), cosmo.dage()
print

cosmo = Cosmology(0.7, 1., 0.)  # h, O_m, O_L
print cosmo.ldist(zvals)
print cosmo.vol_elem(1.), cosmo.dlbt(1.), cosmo.dage()
print

# O_m=0 gives NaN for age:
cosmo = Cosmology(0.7, 0.01, .99)
print cosmo.ldist(zvals)
print cosmo.vol_elem(1.), cosmo.dlbt(1.), cosmo.dage()

zvals = linspace(0.01, 10., 200)

def get_vol(h, Om, Ol, zvals):
    vvals = zeros(len(zvals))
    cosmo = Cosmology(h, Om, Ol)
    for i, z in enumerate(zvals):
        vvals[i] = cosmo.vol_elem(z)
    return vvals

vvals = get_vol(.7, .3, .7, zvals)
loglog(zvals, vvals)
vvals = get_vol(.7, .25, .7, zvals)
loglog(zvals, vvals)
vvals = get_vol(.7, .2, .7, zvals)
loglog(zvals, vvals)
#show()

sc = StdCdlCosmology(0.7, 0.01, .99, 1.)  # h, O_m, O_L, F_fid
sc.lum = 1.

