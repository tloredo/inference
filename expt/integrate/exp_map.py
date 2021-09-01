"""
Run some examples comparing integrating a function with exponential
behavior directly, and mapped to an exponential argument.

Created Apr 2009 by Tom Loredo
"""

from pylab import *
from numpy import *
from ccf import ClenshawCurtis, Fejer1

# Set the inverse scale factor for exponential scaling.
alpha = .5

# Set the range to explore.
range = 30.

def y_x(x, alpha=alpha):
    return exp(alpha*x)/alpha

def x_y(y, alpha=alpha):
    return log(alpha*y)/alpha

def rexp(x, range=range):
    """
    Rolling exponential, -.6 at 1, -.4 at 10
    """
    return exp((alpha + .2*(x-.5*range)/range)*x)

def rexp_map(y, alpha=alpha, range=range):
    """
    Rolling exponential, transformed.
    """
    x = x_y(y)
    return exp((alpha + .2*(x-.5*range)/range)*x) / exp(alpha*x)

l, u = 1., range
y_l, y_u = y_x(l), y_x(u)
for n in [5, 10, 15, 25, 50, 100, 200, 500]:
    ccx = ClenshawCurtis(n, l, u)
    ccy = ClenshawCurtis(n, y_l, y_u)
    xq = ccx.quad(rexp)
    yq = ccy.quad(rexp_map)
    print n, xq, yq

semilogy(ccx.nodes, rexp(ccx.nodes), 'b.')
xlabel('$x$')
ylabel('$f(x)$')
figure()
semilogy(ccy.nodes, rexp_map(ccy.nodes), 'g.')
cc = ClenshawCurtis(25, y_l, y_u)
semilogy(cc.nodes, rexp_map(cc.nodes), 'r.')
xlabel('$y$')
ylabel(r'$f(x)/\exp(\alpha x)$')
show()
