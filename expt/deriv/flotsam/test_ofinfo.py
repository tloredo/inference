from numpy import *
from _ofinfo import pk_scales, ninfol

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
print '*** Testing pk_scales ***'
print 'pk_scales inputs:', mu, maxll, dx
dx, status = pk_scales(llike, mu, maxll, 1., dx, 20, 0.01)
print 'Status should be 0:'
print 'status:', status
print 'index, dx, x, llike (should be 1):'
for i in range(len(dx)):
	x = mu0.copy()
	x[i] += dx[i]
	print i, dx[i], x, llike(x)
print

niter = 5
#ws1, ws2 = zeros(niter,Float), zeros(niter,Float)
info,ierr = ninfol(llike, mu, maxll, dx, niter)

print '*** Testing ninfol, 2x2 ***'
print 'Unit sigmas:'
print info
print ierr
print

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
print 'pk_scales inputs:', mu, maxll, dx
dx, status = pk_scales(llike2, mu, maxll, 1., dx, 20, 0.01)
print 'status:', status
for i in range(len(dx)):
	x = mu0.copy()
	x[i] += dx[i]
	print i, dx[i], x, llike2(x)
print
niter = 1
#ws1, ws2 = zeros(niter,Float), zeros(niter,Float)
info,ierr = ninfol(llike2, mu, maxll, dx, niter)

print 'info, err:'
print info
print ierr
