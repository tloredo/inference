# Obsolete; now incorporated into cub2d.

"""
2-D cubatures on the [-1,1] square, from Ronald Cool's online encyclopedia
of cubature formulas.  See:

R. Cools. An Encyclopaedia of Cubature Formulas J. Complexity, 19: 445-453, 2003.
http://dx.doi.org/10.1016/S0885-064X(03)00011-6
http://www.cs.kuleuven.be/~nines/research/ecf/

These all have positive weights; some rules have nodes strictly inside 
the square, others have some nodes on the boundary.  These correspond to
Cool's "quality" designations PI and PB (positive-interior and 
positive-boundary).

Rule struct indicates the #s of generators with these structures:
  (0,0)  (a,0)  (a,a)  (a,b)

Generator meanings:

Generator: [ Cornerpoints of the unit-cube ]
    ( +- 1, +- 1, ..., +- 1 ).
Generator: [ Intersections of the unit-cube and the axes ]
    ( +- 1, 0, ..., 0 ) and all permutations.
Generator: [ Origin ]
    ( 0, 0, ..., 0 ).
Generator: [ Permutation symmetry ]
    If the cubature formula contains the point ( x1, x2, ... ,xn ) then the
    formula also contains the point ( xp1, xp2, ... ,xpn ) where ( p1, p2, ... ,pn )
    is any permutation of ( 1, 2, ... , n ) and all the coefficients ( weights )
    of each of the points in the cubature formula are the same. If all xi are
    distinct, we have n!points.
Generator: [ Fully symmetric ]
    If the cubature formula contains the point ( x1, x2,... ,xn ) then the formula
    also contains the point ( +- xp1, +- xp2, ... ,+- xpn ) and all the coefficients
    ( weights ) of each of the points in the cubature formula are the same. If all
    xi are distinct and nonzero, we have 2n (n!) points.
Generator: [ Circular symmetric ]
    The nodes and corresponding weights are given in a table. All rows are
    numbered and besides that number there are three more columns. The first
    denotes the number of points that lie on the same circle ( or N-gon ). The
    second column contains the radius of the circle. The weight of those
    points can be found in the last column.
"""

from numpy import *

# Cool's rule c2-3-4a:
gen_3d4pt = ( 0.816496580927726032732428024901963, 0. )

# Cool's rule c2-7-12a:
gen_7d12pt = [
    (array([0.925820099772551461566566776583999, 0.]), 0.241975308641975308641975308641975),
    (array([0.380554433208315656379106359086394, 0.380554433208315656379106359086394]), 0.520592916667394457139919432046731),
    (array([0.805979782918598743707856181350744, 0.805979782918598743707856181350744]), 0.237431774690630234218105259311293)]

def square7d12pt(f):
    """
    Cubature of f(x,y) on the unit square, of degree 7, using 12 points.
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

def square1d4pt(f):
    """
    Cubature of f(x,y) on a rectangle, of degree 1, using 4 points.
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
    Cubature of f(x,y) on a rectangle, of degree 3, using 4 points.
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

if 0:
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

