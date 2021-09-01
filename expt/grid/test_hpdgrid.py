from numpy import *
from inference.grid.hpdgrid import HPDGrid1D, HPDGrid2D
from inference.grid.hpd import HPD1D, HPD2D
import pylab as pl

# A 2-d normal distribution:
norm2 = log(2*pi)


def log_norm2(x,y):
    return -0.5*(x**2 + y**2) - norm2


# A 1-d normal distribution:
norm1 = 0.5*norm2


def log_norm1(x):
    return -0.5*(x**2) - norm1

# A 1-d Gaussian, unnormalize:


def log_gauss1(x):
    return -0.5*(x**2) - norm1


nx = ny = 51

print('1-D normal tests:')
hpd1 = HPDGrid1D(log_norm1, -5., 5., nx)
hpd1_2 = HPD1D(hpd1.logpdf, hpd1.xlo, hpd1.xhi)
print('1-sig:', hpd1.critlevel(.683))
print('      ', hpd1_2.critlevel(.683))
print('2-sig:', hpd1.critlevel(.954))
print('      ', hpd1_2.critlevel(.954))
print()

print('1-D Gaussian tests:')
hpd1g = HPDGrid1D(log_gauss1, -5., 5., nx)
hpd1g_2 = HPD1D(hpd1.logpdf, hpd1.xlo, hpd1.xhi)
print('1-sig:', hpd1g.critlevel(.683))
print('      ', hpd1g_2.critlevel(.683))
print('2-sig:', hpd1g.critlevel(.954))
print('      ', hpd1g_2.critlevel(.954))
print()

print('2-D normal tests:')
hpd2 = HPDGrid2D(log_norm2, -5., 5., nx, -5., 5., ny)
hpd2_2 = HPD2D(hpd2.logpdf, hpd2.xlo, hpd2.xhi, hpd2.ylo, hpd2.yhi)
print('1-sig:', hpd2.critlevel(.683))
print('      ', hpd2_2.critlevel(.683))
print('2-sig:', hpd2.critlevel(.954))
print('      ', hpd2_2.critlevel(.954))
