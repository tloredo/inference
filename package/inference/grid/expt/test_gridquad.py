from scipy import *
from numutils import linEnumStepper
from gridquad import *
from hpd import HPDGrid1d, HPDGrid2d

norm2 = log(2*pi)
def logGauss2(x,y):
	return -0.5*(x**2 + y**2) - norm2

norm1 = 0.5*norm2
def logGauss1(x):
	return -0.5*(x**2) - norm1

nx = ny = 51
dx = dy = 10./(nx-1)

z1 = zeros((nx,), float)
z2 = zeros((nx, ny), float)

for i, x in linEnumStepper(-5., 5., nx):
	z1[i] = logGauss1(x)
	for j, y in linEnumStepper(-5., 5., ny):
		z2[i,j] = logGauss2(x,y)

## z2 = exp(z2)
## print 'To -10:', qgt2d(z2, exp(-norm2-10.))*dx*dy
## print 'To -4:', qgt2d(z2, exp(-norm2-4.))*dx*dy

hpd1 = HPDGrid1d(z1, -5., 5.)
hpd2 = HPDGrid2d(z2, -5., 5., -5., 5.)
