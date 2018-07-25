"""
2-D cubature using monomial rules compiled by Ronald Cools at:

Encyclopaedia of Cubature Formulas (ECF)
http://www.cs.kuleuven.be/~nines/research/ecf/

This module implements selected rules from the ECF:

* Rules on the [-1,1] square with positive weights, some with nodes strictly
  inside the square (Cool's quality designation "PI"), and some that include
  nodes on the boundary ("PB").

* The rules above, mapped to a rectangle.

* Rules for the whole plane with exp(-r^2) weight (built into the rule
  weights).

Rule names (in rule_dict, a dictionary of abscissas and weights) are of the 
following form:

    <region><weight function>_<degree>_<points>

optionally followed by a letter ('a', 'b' etc.) if more than one rule
with the same specifications is offered.  

<region> here is "c2" for the square (2-D cube) and e2 for the whole plane (the
space E_2).

<weight function> here is r2 for exp(-r^2) weighting and absent for no
weighting, e.g, for the c2 rules.  For E_2, the ECF also offers exp(-r),
but these are not yet implemented here.

The rules exactly integrate sums of monomials (x^m * y^n) with degree
m+n up to <degree>, using a total number of function evaluations given
by <points>.  

For the E_2 rules, to give you some idea of the region spanned by the nodes
for the various degrees, here are rough specs for the rules included here:

Degree  Largest abscissa
3       0.8
5       1.6
7       1.7
9       2.3, 2.1
11      2.8, 2.9
13      3.0, 3.1

Selected rules are implemented as functions (besides the abscissa
and weight arrays in rule_dict).

To cite these rules, refer to the web site, which makes the following
request:

When material presented here is used, please refer to the following 
published articles:

R. Cools. 
An Encyclopaedia of Cubature Formulas 
J. Complexity, 19: 445-453, 2003. 
Electronic version:  http://dx.doi.org/10.1016/S0885-064X(03)00011-6

R. Cools. 
Monomial cubature rules since ``Stroud'': a compilation -- part 2. 
J. Comput. Appl. Math., 112(1-2): 21--27, 1999. 
Electronic version:  http://dx.doi.org/10.1016/S0377-0427(99)00229-0

R. Cools and P. Rabinowitz. 
Monomial cubature rules since ``Stroud'': a compilation. 
J. Comput. Appl. Math., 48:309-326, 1993.

This implementation is authored by Tom Loredo and is released subject to 
the SciPy license.
"""

from numpy import pi, sqrt, array

# Import rules auto-generated from ECF web output; currently these are e2 rules.
# Rules on the square are simple and coded directly in this file.
from _cub2d_auto import rule_dict

__all__ = ['rule_dict', 'NormWtCub2D', 'rect3d4pt', 'rect7d12pt', 
    'square1d4pt', 'square3d4pt', 'square7d12pt']


#==============================================================================
# Useful constants
#==============================================================================

rt2 = sqrt(2)

# Cool's rule c2-3-4a abscissas:
gen_3d4pt = ( 0.816496580927726032732428024901963, 0. )

# Cool's rule c2-7-12a abscissas & weights:
gen_7d12pt = [
    (array([0.925820099772551461566566776583999, 0.]), 0.241975308641975308641975308641975),
    (array([0.380554433208315656379106359086394, 0.380554433208315656379106359086394]), 0.520592916667394457139919432046731),
    (array([0.805979782918598743707856181350744, 0.805979782918598743707856181350744]), 0.237431774690630234218105259311293)]


#==============================================================================
# Integrators derived from the "raw" 2-D rules
#==============================================================================

class NormWtCub2D(object):
    """
    2-D cubature using a normal distribution weight function.
    
    The default rule is 'e2r2_7_12' of degree 7 with 12 nodes.  For other
    choices, examine the keys of the dictionary rule_dict in this module.
    See the module documentation for details about the rules.
    """

    def __init__(self, locn, scale, rule='e2r2_7_12'):
        self.rule = rule
        func, self.absc, self.wts = rule_dict[rule]
        self.absc, self.wts = self.absc.copy(), self.wts.copy()  # don't modify originals!
        if len(locn) != 2 or len(scale) != 2:
            raise ValueError('locn and scale must be of length 2!')
        self.absc[:,0] = locn[0] + rt2*scale[0]*self.absc[:,0]
        self.absc[:,1] = locn[1] + rt2*scale[1]*self.absc[:,1]
        self.wts = self.wts / pi

    def __call__(self, func):
        cub = 0.
        for absc, wt in zip(self.absc, self.wts):
            cub += wt * func(*absc)
        return cub


def rect3d4pt(f, xrange, yrange):
    """
    Cubature of f(x,y) on a rectangle, of degree 3, using 4 points.
    All nodes are strictly inside the rectangle.
    """
    xl, xu = xrange
    xm = 0.5*(xu+xl)
    hx = 0.5*(xu-xl)
    yl, yu = yrange
    ym = 0.5*(yu+yl)
    hy = 0.5*(yu-yl)
    px = gen_3d4pt[0]*hx
    py = gen_3d4pt[0]*hy
    cub = f(px+xm,ym) + f(-px+xm,ym) + f(xm,py+ym) + f(xm,-py+ym)
    return cub*hx*hy

def rect7d12pt(f, xrange, yrange):
    """
    Cubature of f(x,y) on a rectangle, of degree 7, using 12 points.
    All nodes are strictly inside the rectangle.
    """
    xl, xu = xrange
    xm = 0.5*(xu+xl)
    hx = 0.5*(xu-xl)
    yl, yu = yrange
    ym = 0.5*(yu+yl)
    hy = 0.5*(yu-yl)
    g, w = gen_7d12pt[0]
    px = g[0]*hx
    py = g[0]*hy
    cub = w * (f(px+xm,ym) + f(-px+xm,ym) + f(xm,py+ym) + f(xm,-py+ym))
    g, w = gen_7d12pt[1]
    px = g[0]*hx
    py = g[0]*hy
    cub += w * (f(px+xm,py+ym) + f(px+xm,-py+ym) + f(-px+xm,py+ym) + f(-px+xm,-py+ym))
    g, w = gen_7d12pt[2]
    px = g[0]*hx
    py = g[0]*hy
    cub += w * (f(px+xm,py+ym) + f(px+xm,-py+ym) + f(-px+xm,py+ym) + f(-px+xm,-py+ym))
    return cub*hx*hy


#==============================================================================
# 2-D monomial rules for the unit [-1,1] square
#==============================================================================

def square1d4pt(f):
    """
    Cubature of f(x,y) on the unit [-1,1] square, of degree 1, using 4 points.
    All nodes are on the corners.
    
    This uses rule c2-1-4 from Ronald Cool's on-line Encyclopedia of
    Cubature Formulas.  This is a trivial rule that simply sums the
    corner contributions with unit weights.

    Region: Cube 
    Dimension: 2 
    Degree: 1 
    Points: 4 
    Structure: Fully symmetric 
    Rule struct: 0 0 1 0 
    Product trapezoidal formula
    Generator: [ Cornerpoints of the unit-cube ] 
    ( 1., 1., ) 
    Corresponding weight: 
    1.,
    """
    return f(1.,1.) + f(1.,-1.) + f(-1.,1.) + f(-1.,-1.)


def square3d4pt(f):
    """
    Cubature of f(x,y) on the unit [-1,1] square, of degree 3, using 4 points.
    All nodes are strictly inside the square.

    This uses rule c2-3-4a from Ronald Cool's on-line Encyclopedia of
    Cubature Formulas.

    Region: Cube 
    Dimension: 2 
    Degree: 3 
    Points: 4 
    Structure: Fully symmetric 
    Rule struct: 0 1 0 0 
    Generator: [ Fully symmetric ] 
    ( 0.816496580927726032732428024901963, 0., ) 
    Corresponding weight: 
    1.,
    """
    p = gen_3d4pt[0]
    return f(p,0.) + f(-p,0.) + f(0.,p) + f(0.,-p)

def square7d12pt(f):
    """
    Cubature of f(x,y) on the unit [-1,1] square, of degree 7, using 12 points.
    All nodes are strictly inside the square.

    This uses rule c2-7-12a from Ronald Cool's on-line Encyclopedia of
    Cubature Formulas.
    
    Region: Cube 
    Dimension: 2 
    Degree: 7 
    Points: 12 
    Structure: Fully symmetric 
    Rule struct: 0 1 2 0 
    The points are the common zeros of the orthogonal polynomials
    P4,1 = x3y - xy3
    
    P4,2 = x4 - y4 - 6/7x2 + 6/7y2
    
    P4,3 = x4 + 54/55x2y2 + y4 - 456/385x2 - 456/385y2 + 108/385
    
    
    Generator: [ Fully symmetric ] 
    ( 0.925820099772551461566566776583999, 0., ) 
    Corresponding weight: 
    0.241975308641975308641975308641975,
    Generator: [ Fully symmetric ] 
    ( 0.380554433208315656379106359086394, 0.380554433208315656379106359086394, ) 
    Corresponding weight: 
    0.520592916667394457139919432046731,
    
    Generator: [ Fully symmetric ] 
    ( 0.805979782918598743707856181350744, 0.805979782918598743707856181350744, ) 
    Corresponding weight: 
    0.237431774690630234218105259311293,
    """
    g, w = gen_7d12pt[0]
    p = g[0]
    cub = w * (f(p,0.) + f(-p,0.) + f(0.,p) + f(0.,-p))
    g, w = gen_7d12pt[1]
    p = g[0]
    cub += w * (f(p,p) + f(p,-p) + f(-p,p) + f(-p,-p))
    g, w = gen_7d12pt[2]
    p = g[0]
    cub += w * (f(p,p) + f(p,-p) + f(-p,p) + f(-p,-p))
    return cub


#==============================================================================
# Some yet-to-be implemented rule descriptions
#==============================================================================

"""
Rule c2-3-5, quality PI
Region: Cube 
Dimension: 2 
Degree: 3 
Points: 5 
Structure: Fully symmetric 
Rule struct: 1 1 0 0 
Generator: [ Intersections of the unit-cube and the axes ] 
( 1., 0., ) 
Corresponding weight: 
0.666666666666666666666666666666666,
Generator: [ Origin ] 
( 0., 0., ) 
Corresponding weight: 
1.33333333333333333333333333333333,
"""

"""
Rule c2-3-5a, quality PI
Region: Cube 
Dimension: 2 
Degree: 3 
Points: 5 
Structure: Fully symmetric 
Rule struct: 1 0 1 0 
Generator: [ Cornerpoints of the unit-cube ] 
( 1., 1., ) 
Corresponding weight: 
0.333333333333333333333333333333333,
Generator: [ Origin ] 
( 0., 0., ) 
Corresponding weight: 
2.66666666666666666666666666666666,
"""
