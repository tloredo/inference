from numpy import *
from inference.integrate.cub2d import rule_dict, NormWtCub2D

def unity(x,y):
    """
    Integral of unity with exp(-r^2) weight should be pi.
    """
    return 1.

def xval(x,y):
    """
    Integral of xval with exp(-r^2) weight should be zero.
    """
    return x

def yval(x,y):
    """
    Integral of yval with exp(-r^2) weight should be zero.
    """
    return y

def xsqr(x,y):
    """
    Integral of xsqr with exp(-r^2) weight should be pi/2.
    """
    return x*x

def ysqr(x,y):
    """
    Integral of ysqr with exp(-r^2) weight should be pi/2.
    """
    return y*y


locn, scale = array([200., 0.01]), array([2., 20.5])

def test_rule(name):
    cub = NormWtCub2D(locn, scale, rule=name)
    
    norm = cub(unity)
    xbar, ybar = cub(xval), cub(yval)
    print 'Testing with', name
    print '  Norm: ', norm, '(', norm-1., ')'
    print '  Means: ', xbar, ybar, '(', xbar-locn[0], ybar-locn[1], ')'
    x2, y2 = cub(xsqr), cub(ysqr)
    xsd, ysd = sqrt(x2-xbar**2), sqrt(y2-ybar**2)
    print '  Std deviations: ', xsd, ysd, '(', xsd-scale[0], ysd-scale[1], ')'
    print

for name in rule_dict.keys():
    test_rule(name)


# Test square rules:

if 1:
    from inference.integrate.cub2d import square7d12pt, rect7d12pt, square3d4pt, rect3d4pt
    print 'Testing constant integrand:'
    def g(x,y):
        return 1
    print square7d12pt(g), '(should be 4.0)'
    print rect7d12pt(g, (0,1), (0,1)), '(should be 1.0)'

    print 'Testing 2nd degree integrand:'
    def f(x,y):
        # print x,y
        return 1 + x + x*y + 2*y
    print square7d12pt(f), '(should be 4.0)'
    print rect7d12pt(f, (0,1), (0,1)), '(should be 2.75)'
    print rect3d4pt(f, (0,1), (0,1)), '(should be 2.75)'

    print 'Testing 2nd degree integrand:'
    def f(x,y):
        # print x,y
        return 1 + x**2 + y**2
    print square7d12pt(f), '(should be ?)'
    print square3d4pt(f), '(should be ?)'

    print 'Testing 3rd degree integrand:'
    def f(x,y):
        # print x,y
        return 1 + 3*y**2 + 0.5*y*x**2 + x**3
    print square7d12pt(f), '(should be 8)'
    print square3d4pt(f), '(should be 8)'
    print rect7d12pt(f, (0,1), (0,1)), '(should be 2 1/3)'
    print rect3d4pt(f, (0,1), (0,1)), '(should be 2 1/3)'

    print 'Testing 4th degree integrand:'
    def f(x,y):
        # print x,y
        return 1 + x + 3*(x*y)**2
    print square7d12pt(f), '(should be 5 1/3)'
    print square3d4pt(f), '(should differ from 5 1/3)'
    print rect7d12pt(f, (0,1), (0,1)), '(should be 1 5/6 = 1.8333...)'
    print rect3d4pt(f, (0,1), (0,1)), '(should differ from 1 5/6)'
