"""
Run some examples comparing integrating a function with power law
behavior directly, and mapped to a power law argument.

Created Apr 2009 by Tom Loredo
"""

from pylab import *
from numpy import *
from ccf import ClenshawCurtis, Fejer1

# Set the power law for the map; x**(-beta).
beta = 2.5
beta1 = 1 - beta
assert beta != 1.

# Set the range to explore.
range = 100.
mid = sqrt(range)

def y_x(x, beta1=beta1):
    return x**beta1 / beta1

def x_y(y, beta1=beta1):
    logx = log(beta1*y)/beta1
    return exp(logx)

def roll(x, beta=beta, mid=mid):
    """
    Rolling power law.
    """
    pl1 = -beta + .5
    return (x/mid)**pl1 / (1. + x/mid)

def roll_map(y, beta=beta, range=range, mid=mid):
    """
    Rolling exponential, transformed.
    """
    x = x_y(y)
    return roll(x) * x**beta

l, u = 1., range
y_l, y_u = y_x(l), y_x(u)
for n in [5, 10, 15, 25, 50, 100, 200, 500]:
    ccx = ClenshawCurtis(n, l, u)
    ccy = ClenshawCurtis(n, y_l, y_u)
    xq = ccx.quad(roll)
    yq = ccy.quad(roll_map)
    print n, xq, yq

loglog(ccx.nodes, roll(ccx.nodes), 'b.')
cc = ClenshawCurtis(25, l, u)
loglog(cc.nodes, roll(cc.nodes), 'r.')
xlabel('$x$')
ylabel('$f(x)$')
figure()
loglog(-ccy.nodes, roll_map(ccy.nodes), 'g.')
cc = ClenshawCurtis(25, y_l, y_u)
loglog(-cc.nodes, roll_map(cc.nodes), 'r.')
xlabel('$-y$')
ylabel(r'$f(x)/x^{-\beta}$')
show()
