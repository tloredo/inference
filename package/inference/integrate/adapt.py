import _adapt

"""
For adaptmodule.adapt:

adapt - Function signature:
  minpts,relerr,finest,ifail = adapt(a,b,minpts,maxpts,functn,eps,lenwrk,[ndim,functn_extra_args])
Required arguments:
  a : input rank-1 array('d') with bounds (ndim)
  b : input rank-1 array('d') with bounds (ndim)
  minpts : input int
  maxpts : input int
  functn : call-back function => func
  eps : input float
  lenwrk : input int
Optional arguments:
  ndim := len(a) input int
  functn_extra_args := () input tuple
Return objects:
  minpts : int
  relerr : float
  finest : float
  ifail : int
Call-back functions:
  def func(args,[n]): return r
  Required arguments:
    args : input rank-1 array('d') with bounds (n)
  Optional arguments:
    n := len(args) input int
  Return objects:
    r : float
"""

def adapt_sizes(ndim, maxfac=100):
    """
    Integers constraining the adapt algorithm: minpts, maxpts, lenwrk.

    From ADAPT's documentation:

    MAXPTS  MAXIMUM NUMBER OF FUNCTION EVALUATIONS TO BE ALLOWED,
            WHICH MUST BE AT LEAST RULCLS, WHERE
            RULCLS =  2**NDIM+2*NDIM**2+6*NDIM+1

              FOR NDIM =  2   3   4   5   6   7   8   9   10
              MAXPTS >=  25  45  73 113 173 269 433 729 1285
           A suggested value for MAXPTS is 100 times the above values.

    LENWRK  LENGTH OF ARRAY WRKSTR OF WORKING STORAGE, THE ROUTINE
            NEEDS (2*NDIM+3)*(1+MAXPTS/RULCLS)/2 FOR LENWRK IF
            MAXPTS FUNCTION CALLS ARE USED.
            FOR GUIDANCE, IF YOU SET MAXPTS TO 100*RULCLS (SEE TABLE
            ABOVE) THEN ACCEPTABLE VALUES FOR LENWRK ARE

              FOR NDIM = 2    3    4    5    6    7    8     9
              LENWRK =  357  561  1785 3417 6681 13209 26265 52377

    Note:  The LENWRK expression significantly underestimates the
    LENWRK values given above.
    """
    rulcls = 2**ndim + 2*ndim**2 + 6*ndim + 1
    minpts, maxpts = rulcls, maxfac*rulcls
    lenwrk = (2*ndim + 3) * (1 + maxfac) / 2
    return minpts, maxpts, lenwrk

def adapt(func, lower, upper, relacc, maxfac=100, wrkfac=1):
    """
    Adaptive quadrature over a hyper-rectangle using Genz's ADAPT algorithm.

    func = function to be integrated, func(x), for x a 1-d Float array
    lower, upper = arrays with lower and upper bounds on x
    relacc = required relative accuracy
    maxfac = factor by which the max # of points used is allowed to
             exceed the min #
    wrkfac = factor controlling the amount of working space reserved
             for the calculation

    Returns:  result, err, npts.

      result = estimated integral
      err = estimate of relative error of result
      ntps = # of function evaluations used
    """
    minpts, maxpts, lenwrk = adapt_sizes(len(lower))
    npts, err, result, ifail = _adapt.adapt(lower, upper, minpts, maxpts, \
                                     func, relacc, lenwrk)
    if ifail==0:
        return result, err, npts
    elif ifail==1:
        raise RuntimeError, 'Failed to achieve specified accuracy; increase maxfac!'
    elif ifail==2:
        raise RuntimeError, 'Failed to achieve specified accuracy; increase wrkfac!'
