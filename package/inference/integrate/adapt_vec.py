from numpy import zeros
from ._dcuhre import dcuhre

__all__ = ['adapt_vec']


def adapt_vec(nfuncs, vfunc, lower, upper, releps,
              abseps=0., key=0, minpts=0, maxpts=100000):
    """Adaptive cubature of a vector of functions over a hyperrectangle.

    Successful completion returns (result, abserr, neval):
        result  =   vector of integration results
        abserr  =   vector of absolute errors for integrals
        neval   =   number of function evaluations used
    A successful return abserr <=  abseps or abserr <=  abs(result)*releps 
    with maxpts or fewer function evaluations, for all integrals.

    The DCUHRE subregion-adaptive cubature algorithm of Bernsten, Espelid & Genz.
    Double-precision Adaptive CUbature over Hyper REctangles

    J.Berntsen, T.O.Espelid and A.Genz, An Adaptive Algorithm
    for the Approximate Calculation of Multiple Integrals,

    J.Berntsen, T.O.Espelid and A.Genz, DCUHRE: An Adaptive
    Multidimensional Integration Routine for a Vector of Integrals
    """

    ndim = len(lower)
    if ndim < 2:
        raise ValueError('Dimension must be >=2!')
    elif ndim > 15:
        raise ValueError('Dimension limited to [2,15]!')

    # Calculate the number of vfunc evaluations per region.
    # Min allowed maxpts is 3*num.
    if (key == 0 or key == 1) and ndim == 2:
        num = 65
    elif (key == 0 or key == 2) and ndim == 3:
        num = 127
    elif (key == 0 and ndim >= 3) or key == 3:
        num = 1 + 8*ndim + 2*ndim*(ndim-1) + 4*ndim*(ndim-1) + \
            4*ndim*(ndim-1)*(ndim-2)/3 + 2**ndim
    elif key == 4:
        num = 1 + 6*ndim + 2*ndim*(ndim-1) + 2**ndim
    else:
        raise ValueError('Illegal key value!')

    # Calculate the max number of subregions, and workspace size.
    maxsub = (maxpts-num)/(2*num) + 1
    nw = maxsub*(2*ndim+2*nfuncs+2) + 17*nfuncs + 1
    work = zeros(nw, float)
    print('Rule size, max subregions, workspace size:', num, maxsub, nw)

    # Store result in a persistant array, potentially allowing restarts.
    result = zeros(nfuncs, float)
    restart = 0

    result, abserr, neval, ifail = dcuhre(nfuncs, lower, upper, minpts,
                                          maxpts, vfunc, abseps, releps, key, restart, result, work)

    if ifail == 0:
        return result, abserr, neval
    elif ifail == 1:
        print('Inaccurate result:', result, abserr, neval)
        raise ValueError('maxpts too small to achieve requested accuracy!')
    elif ifail == 2:
        raise ValueError('Illegal key value!')
    elif ifail == 3:
        raise ValueError('Illegal dimension!')
    elif ifail == 4 or ifail == 5:
        raise ValueError('Illegal key for integrand dimension!')
    elif ifail == 6:
        raise ValueError('Number of integrands == 0!')
    elif ifail == 7:
        raise ValueError('Input region has zero volume!')
    elif ifail == 8:
        raise ValueError('maxpts less than 3 * rule size!')
    elif ifail == 9:
        raise ValueError('maxpts is < minpts!')
    elif ifail == 10:
        raise ValueError('Requested accuracy < 0!')
    elif ifail == 11:
        raise ValueError('Workspace too small!')
    elif ifail == 12:
        raise ValueError('Illegal restart value!')
    else:
        raise RuntimeError('Unknown exit code from dcuhre!')
