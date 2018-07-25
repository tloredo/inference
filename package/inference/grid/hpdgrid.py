"""
Objects for calculating highest posterior density credible regions.

HPDGrid1D and HPDGrid2D find the boundaries of HPD regions by evaluating a
(log) probability density function over a 1-D or 2-D grid, and then integrating
over the grid.  The function gets (re)normalized over the grid; e.g., one can
use an unnormalized prior*likelihood as the pdf (the normalization constant
is accessible as a data attribute and may be useful for Bayesian model
comparison).

For situations where other tools have previously evaluated a pdf over a grid,
see the classes of the same name in the hpd module.

By Tom Loredo
"""

from inference.utils.numutils import zbracket, lin_stepper, log_stepper
from scipy import exp, log, sqrt, pi, array, stats, zeros
from scipy.optimize import brentq
from ._gridquad import qgt1d, qgt2d, xvalues
import pylab as pl

__all__ = ['HPDGrid1D', 'HPDGrid2D']

rt2pi = sqrt(2*pi)
stdnorm = stats.norm()


def sig2prob(nd, nsig):
    """
    The probability within nsig std deviations of the mean of a n-dim normal
    distribution.
    """
    # For 1d the upper tail area is given by the std normal survival function.
    #utail = stdnorm.sf(nsig)
    # return 1. - 2*utail
    return stats.chi2.cdf(nsig**2, nd)


def prob2sig(nd, p):
    """
    The number of std deviations corresponding to probability p for n-dim normal.
    """
    # These two lines are for the 1-d normal case.
    #utail = 0.5*(1.-p)
    # return stdnorm.isf(utail)
    return sqrt(stats.chi2.ppf(p, nd))


def prob2dens(nd, p):
    """
    The n-dim normal PDF density bounding an HPD credible region of content p.
    """
    nsig = prob2sig(nd, p)
    return exp(-0.5*nd * nsig**2)/rt2pi**nd


def prob2ratio(nd, p):
    """
    The relative normal PDF density bounding an HPD credible region
    of content p (relative to the peak density).
    """
    nsig = prob2sig(nd, p)
    return exp(-0.5*nd * nsig**2)


def sig2ratio(nd, nsig):
    """
    The relative normal PDF density bounding an HPD credible region
    spanning nsig std deviations (relative to the peak density).
    """
    return exp(-0.5*nd * nsig**2)


class HPDGrid1D(object):
    """
    Evaluate a 1-D PDF on a grid and calculate HPD region boundaries.

    For a pdf p(x), this object calculates and stores p(x) values across
    of x values spanning [xlo, xhi] that are either linearly spaced, or
    linearly spaced in log(x).  HPD regions may subsequently be calculated
    using the grid.  For grids linear in x, the regions are for p(x).  For
    grids linear in log(x), the regions are for p(log x), but with boundaries
    specified in terms of x itself, appropriate for use on plots of 
    p(log x) (equal to x*p(x)) vs. x with a logarithmic abscissa.
    """

    # Probabilities associated with normal dist'n "sigmas":
    p1sig, p2sig, p3sig = sig2prob(1,1), sig2prob(1,2), sig2prob(1,3)

    # Corresponding relative PDF levels defining normal HPD regions,
    # to use as starting points:
    ratio1 = sig2ratio(1,1)
    ratio2 = sig2ratio(1,2)
    ratio3 = sig2ratio(1,3)

    def __init__(self, log_pdf, xlo, xhi, nx, step='linear'):
        """
        Initialize the HPDGrid1D object.

        :Parameters:

          log_pdf : float function of scalar argument
            Function returning log pdf value log_pdf(x)

          xlo, xhi : float
            Range spanned by logpdf

          nx : int
            Number of points in the logpdf grid

          step : 'linear' OR 'lin' OR 'logarithmic' OR 'log'
            'linear' for a grid linear in x; 'log' for a grid linear in log(x)
        """
        self.log_pdf = log_pdf
        xlo, xhi = float(xlo), float(xhi)
        self.xlo, self.xhi = xlo, xhi
        self.nx = nx
        if step == 'linear' or step == 'lin':
            self.linstep = True
            self.xvals = array([x for x in lin_stepper(xlo, xhi, self.nx)])
            self.dx = (xhi-xlo)/(self.nx-1)
        elif step == 'log' or step == 'logarithmic':
            self.linstep = False
            self.xvals = [x for x in log_stepper(xlo, xhi, self.nx)]
            self.dx = log(xhi/xlo)/(self.nx-1)
        else:
            raise ValueError('Invalid step type!')
        logpdf = zeros(nx, float)
        for i, x in enumerate(self.xvals):
            logpdf[i] = log_pdf(x)
        self.logpdf = logpdf
        self.max = logpdf.max()
        self.spdf = exp(logpdf-self.max)  # scaled PDF
        # For log steps, change variables to log(x).
        if not self.linstep:
            for i, x in enumerate(self.xvals):
                self.spdf[i] = self.spdf[i]*x
        # Find the overall normalization for spdf.
        self.snorm = qgt1d(self.spdf, self.spdf.min())*self.dx
        # lml is the log marginal likelihood if logpdf = prior*like.
        self.lml = log(self.snorm) + self.max
        self.probs = []
        self.deltas = []
        self.levels = []
        self.bounds = []

    def norm(self):
        """
        Return the log of the normalization constant for the grid.

        If the log_pdf function is prior*likelihood, this is the marginal
        likelihood.
        """
        return log(self.snorm) + self.max

    def fracgt(self, logratio):
        """
        The fraction of the posterior with log density > maximum + logratio.
        """
        ratio = exp(logratio)
        return qgt1d(self.spdf, ratio)*self.dx/self.snorm

    def critlevel(self, p=None, snsig=None, tol=3.e-4):
        """
        The log critical PDF ratio (wrt mode) bounding an HPD credible region
        with probability p.
        """
        # Take as an initial guess the ratio for a 1-d normal.
        if snsig:
            if p:
                raise ValueError('Provide only one of p, snsig!')
            lr1 = log(sig2ratio(1, snsig))
            p = sig2prob(1, snsig)
        else:
            lr1 = log(prob2ratio(1, p))
        # p = float(p)  # Old workaround for a now-fixed scalar array bug.
        # Use the guess to bracket the solution and then solve.
        lr2 = lr1 - 0.1
        # print 'Guesses:', lr1, lr2
        self._target = p
        lr1, lr2 = zbracket(self._diff, lr1, lr2)
        # print 'Bracket:', lr1, lr2
        logratio = float(brentq(self._diff, lr1, lr2, xtol=tol))
        # print repr(logratio)
        self.probs.append(p)
        self.deltas.append(logratio)
        self.levels.append(self.max+logratio)
        # 10 is the max # boundaries; it should handle pretty bumpy PDFs!
        bounds, nb, ok = xvalues(self.xvals, self.logpdf-self.max, logratio, 10)
        if not ok:
            raise RuntimeError('Too many boundary points in PDF!')
        bounds = bounds[:nb]
        self.bounds.append(bounds)
        return logratio, bounds

    def _diff(self, logratio):
        """
        The function whose zero gives critRatio; it is fracgt - target prob.
        """
        return self.fracgt(logratio) - self._target

    def plot(self):
        pl.plot(self.xvals, exp(self.logpdf))
        xl, xu = self.xvals[0], self.xvals[-1]
        for level in self.levels:
            plevel = exp(level)
            pl.plot([xl,xu], [plevel,plevel])


class HPDGrid2D(object):

    # Probabilities associated with normal dist'n "sigmas":
    p1sig, p2sig, p3sig = sig2prob(2,1), sig2prob(2,2), sig2prob(2,3)

    # Corresponding relative PDF levels defining normal HPD regions,
    # to use as starting points:
    ratio1 = sig2ratio(2,1)
    ratio2 = sig2ratio(2,2)
    ratio3 = sig2ratio(2,3)

    def __init__(self, log_pdf, xlo, xhi, nx, ylo, yhi, ny, xstep='linear',
                 ystep='linear'):
        self.log_pdf = log_pdf
        self.xlo, self.xhi = xlo, xhi
        self.nx = nx
        if xstep == 'linear' or xstep == 'lin':
            self.xlinstep = True
            self.xvals = array([x for x in lin_stepper(xlo, xhi, self.nx)])
            self.dx = (xhi-xlo)/(self.nx-1)
        elif xstep == 'log' or xstep == 'logarithmic':
            self.xlinstep = False
            self.xvals = array([x for x in log_stepper(xlo, xhi, self.nx)])
            self.dx = log(xhi/xlo)/(self.nx-1)
        else:
            raise ValueError('Invalid xstep type!')
        self.ylo, self.yhi = ylo, yhi
        self.ny = ny
        if ystep == 'linear' or ystep == 'lin':
            self.ylinstep = True
            self.yvals = array([y for y in lin_stepper(ylo, yhi, self.ny)])
            self.dy = (yhi-ylo)/(self.ny-1)
        elif ystep == 'log' or ystep == 'logarithmic':
            self.ylinstep = False
            self.yvals = array([y for y in log_stepper(ylo, yhi, self.ny)])
            self.dy = log(yhi/ylo)/(self.ny-1)
        else:
            raise ValueError('Invalid ystep type!')
        logpdf = zeros((nx, ny), float)
        for i,x in enumerate(self.xvals):
            for j,y in enumerate(self.yvals):
                logpdf[i,j] = log_pdf(x,y)
        self.logpdf = logpdf
        self.max = logpdf.max()
        self.spdf = exp(logpdf-self.max)  # Scaled PDF
        # For log steps, change variables to log(x or y).
        if not self.xlinstep:
            for i, x in enumerate(self.xvals):
                self.spdf[i,:] = self.spdf[i,:]*x
        if not self.ylinstep:
            for j, y in enumerate(self.yvals):
                self.spdf[:,j] = self.spdf[:,j]*y
        # Find the overall normalization for spdf.
        self.norm = qgt2d(self.spdf, self.spdf.min())*self.dx*self.dy
        # lml is the log marginal likelihood if logpdf = prior*like.
        self.lml = log(self.norm) + self.max
        self.probs = []
        self.deltas = []
        self.levels = []

    def fracgt(self, logratio):
        """
        The fraction of the posterior with log density > maximum + logratio.
        """
        ratio = exp(logratio)
        return qgt2d(self.spdf, ratio)*self.dx*self.dy/self.norm

    def critlevel(self, p=None, snsig=None, tol=3.e-4):
        """
        The log critical PDF ratio (wrt mode) bounding an HPD credible region
        with probability p.
        """
        # Take as an initial guess the ratio for a 1-d normal.
        if snsig:
            if p:
                raise ValueError('Provide only one of p, snsig!')
            lr1 = log(sig2ratio(1, snsig))
            p = sig2prob(1, snsig)
        else:
            lr1 = log(prob2ratio(1, p))
        p = float(p)  # Otherwise it may be a scalar array.
        # Use this to bracket the solution and then solve.
        lr2 = lr1 - 0.1
        # print 'Guesses:', lr1, lr2
        self._target = p
        lr1, lr2 = zbracket(self._diff, lr1, lr2)
        # print 'Bracket:', lr1, lr2
        logratio = float(brentq(self._diff, lr1, lr2, xtol=tol))
        # print repr(logratio)
        self.probs.append(p)
        self.deltas.append(logratio)
        self.levels.append(self.max+logratio)
        return logratio

    def _diff(self, logratio):
        """The function whose zero gives critRatio; it is fracgt - target prob."""
        return self.fracgt(logratio) - self._target

    def plot(self):
        x,y = pl.meshgrid(self.xvals,self.yvals)
        pl.figure()
        pl.contour(x,y,self.logpdf,self.levels)
