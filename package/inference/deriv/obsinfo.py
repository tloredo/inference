"""
obsinfo - Tools for calculating the observed Fisher information matrix
from log-likelihood functions.
"""

from _obsinfo import pk_scales, ninfol

__all__ = ['dllsteps', 'obsinfo']

def dllsteps(loglike, mle, guess, dll=1., maxll=None, maxiter=20, tol=0.01):
    """
    Calculate step sizes for parameters leading to a specified change in
    log-likelihood.

    deltas = dllsteps(loglike, mle, guess, dll=1., maxll=None, niter=10, tol=0.01)

    loglike =   log likelihood function of a 1-d array argument
    mle     =   1-d array of maximum likelihood parameter estimates
    guess   =   a 1-d array of initial guesses for the scales
    dll     =   target change in log likelihood defining scales
    maxll   =   loglike(mle); will be evaluated if maxll=None
    maxiter =   maximum number of iterations
    tol     =   tolerance (% of dll) for convergence
    """
    if maxll is None:
        maxll = loglike(mle)
    deltas, status = pk_scales(loglike, mle, maxll, dll, guess, maxiter, tol)
    if status == 3:
        raise RuntimeError('Failed to converge in solving for dll scale!')
    elif status == 2:
        raise ValueError('Input mle is not at maximum!')
    return deltas


def obsinfo(loglike, mle, deltas, maxll=None, niter=2):
    """
    Calculate the observed Fisher information matrix from a log-likelihood
    function, using Richardson extrapolation to improve the estimate.

    info, err = obsinfo(loglike, mle, steps, maxll=None, niter=2)

    The information matrix is - (d/dx_i)(d/dx_j) log(likelihood), i.e.,
    the negative Hessian, divided by the likelihood value at the mle.

    info    =   information matrix
    err     =   array of estimated errors in info

    Note that loglike is not exponentiated anywhere in this algorithm,
    so it is safe to use functions with very large or small logarithms.

    loglike =   log likelihood function of a 1-d array argument
    mle     =   1-d array of maximum likelihood parameter estimates
    deltas  =   1-d array of initial step sizes for the parameters;
                can be calculated using dllscales().
    maxll   =   loglike(mle); will be evaluated if maxll=None
    niter   =   number of Richardson extrapolation steps to use to
                refine the derivative calculations; if niter=0, derivatives
                from simple 2nd differences are reported, and 
                err = an array of zeros
    
    For each iteration, the step sizes are shrunk by 1/sqrt(2).
    """
    # Note ninfol interprets # of iterations as # of sets of 2nd differences
    # to use, not # of Richardson extrapolation steps (i.e., #=1 corresponds
    # to plain 2nd differences, with no extrapolation).
    if maxll is None:
        maxll = loglike(mle)
    return ninfol(loglike, mle, maxll, deltas, niter+1)
