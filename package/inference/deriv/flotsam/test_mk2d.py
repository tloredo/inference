from scipy import *
from mk2d import mk2dc, mk2df, ninfol

a = mk2dc(4)
b = mk2df(4)

print 'intent(c,out)'
print a
print
print 'intent(out)'
print b
print

mu = array([1., 3.])
def workspace(a):
	"""Return a zero array of the same shape and type as a,
	to use as workspace."""
	return zeros(a.shape, a.dtype)

ws1, ws2 = workspace(mu), workspace(mu)
maxll = 0.
dx = .1*mu
info,ierr,ws1,ws2 = ninfol(mu, maxll, dx, 5, ws1, ws2)
