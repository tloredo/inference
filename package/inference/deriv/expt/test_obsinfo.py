from numpy import *
from inference.deriv import dllsteps, obsinfo

mu = array([1., 3.])
mu0 = mu.copy()

# Note that mu (2nd argument in pk_scales) is modified within
# that function; thus it cannot be used as the (global) true
# mean in llike.

def llike(xy):
	val = -0.5*sum( (xy-mu0)**2 )
	# print xy, val
	return val

dx = array([0.2, 0.4])
maxll = llike(mu)
print '*** Testing dllsteps ***'
print 'dllsteps inputs:', mu, maxll, dx
dx = dllsteps(llike, mu, dx, 1., maxll)
print 'index, dx, x, llike (should be -1):'
for i in range(len(dx)):
	x = mu0.copy()
	x[i] += dx[i]
	print i, dx[i], x, llike(x)
print

niter = 5
#ws1, ws2 = zeros(niter,Float), zeros(niter,Float)
info,ierr = obsinfo(llike, mu, dx, maxll)

print '*** Testing obsinfo, 2x2 ***'
print 'Unit sigmas, no correlation -> info, err:'
print info
print ierr
print

print
print '*** Correlated 2x2 tests ***'
sigmas = array([.3, 5.])
cross = 1.5
print 'Sigmas:', sigmas
print 'squared inverses:', 1/sigmas**2
print 'cross coef:', cross
def llike2(xy):
	deltas = xy-mu0
	val = -0.5*(sum((deltas/sigmas)**2) + 2*cross*prod(deltas))
	return val

dx = array([0.2, 0.4])
maxll = llike2(mu)
print 'dllsteps inputs:', mu, maxll, dx
dx = dllsteps(llike2, mu, dx, 1., maxll)
print 'index, dx, x, llike (should be -1):'
for i in range(len(dx)):
	x = mu0.copy()
	x[i] += dx[i]
	print i, dx[i], x, llike2(x)
print
niter = 1
#ws1, ws2 = zeros(niter,Float), zeros(niter,Float)
info,ierr = obsinfo(llike2, mu, dx, maxll)

print 'info, err:'
print info
print ierr
