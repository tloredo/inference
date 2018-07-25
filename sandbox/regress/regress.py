from numpy import array, sqrt, ones, Float
import _regress


def bistats(x, y):
    """
    Basic statistics summarizing a bivariate distribution.
    
    xbar, xstd, ybar, ystd, covar, rho = bistats(x,y)
    
    xbar, ybar are the sample means
    xstd, ystd are the sample standard deviations (with 1/n
        rather than 1/(n-1) factors in the variance sum)
    covar is the covariance
    rho is Pearson's linear correlation coefficient
    """
    if len(x) != len(y): raise ValueError, 'x, y length mismatch!'
    if len(x) < 2: raise ValueError, 'Need at least n=2!'
    xbar, xstd, ybar, ystd, covar, rho, ierr = _regress.datstt(x,y)
    if ierr == 0:
        raise ValueError, 'Input data has zero Sum(x-xbar)(y-ybar)!'
    xstd = sqrt(xstd)
    ystd = sqrt(ystd)
    return xbar, xstd, ybar, ystd, covar, rho


class BiStats(object):
    """
    Basic statistics summarizing a bivariate distribution.
    
    b = BiStats(x, y, nocopy=None) creates an object with these attributes:
    
    xbar, ybar  =   sample means
    xstd, ystd  =   sample standard deviations (with 1/n rather than 1/(n-1) 
                    factors in the variance sum)
    covar       =   covariance
    rho         =   Pearson's linear correlation coefficient
    
    The input data is copied to x, y attributes, unless nocopy is True,
    in which case the x, y attributes are references to the original data.
    """

    def __init__(self, x, y, nocopy=None):
        if len(x) != len(y): raise ValueError, 'x, y length mismatch!'
        if len(x) < 2: raise ValueError, 'Need at least n=2!'
        # Handle copying first so the datstt call doesn't duplicate
        # any needed conversion of a sequence to an array.
        if nocopy:
            self.x = x
            self.y = y
        else:
            try:
                self.x = x.copy()
                self.y = y.copy()
            except AttributeError:
                self.x, self.y = array(x), array(y)
        self.xbar, self.xstd, self.ybar, self.ystd, \
            self.covar, self.rho, error = _regress.datstt(self.x, self.y)
        if error == 0:
            raise ValueError, 'Input data has zero Sum(x-xbar)(y-ybar)!'
        self.xstd = sqrt(self.xstd)
        self.ystd = sqrt(self.ystd)


class LinRegResult(object):
    """
    Results of a bivariate linear regression calculation.
    
    Basic attributes:
        method  =   String identifying regression method
        slope   =   Slope of regression line
        inter   =   Intercept of regression line
        serr    =   Standard error for slope
        ierr    =   Standard error for intercept
        x, y    =   Copies of data vectors used for calculation
        xerr, yerr    =   Copies of x, y errors used for calculation 
                          (if passed)
    """

    def __init__(self, method, slope, inter, serr, ierr, x, y, xerr=None, yerr=None):
        self.method = method
        self.slope = slope
        self.inter = inter
        self.serr = serr
        self.ierr = ierr
        self.x, self.y = x, y
        self.npts = len(self.x)
        if xerr is None:
            self.xmin, self.xmax = min(x), max(x)
        else:
            self.xmin, self.xmax = min(x-xerr), max(x+xerr)
        if yerr is None:
            self.ymin, self.ymax = min(y), max(y)
        else:
            self.ymin, self.ymax = min(y-yerr), max(y+yerr)

        
    
def wlss(x, y, yerr):
    """
    Perform weighted least squares with scatter for fitting a line to
    y vs. x data, where the y measurements have known measurement
    errors (yerr) but also additional (homoscedastic) intrinsic scatter.
    
    yerr may be a scalar for constant (homoscedastic) measurement error,
    or an array for (heteroscedastic) errors that vary from datum to datum.
    
    Returns a LinRegResult object with an additional attribute:
        scat    =   Standard deviation of the intrinsic scatter
        yerr    =   Copy of y errors used for calculation
    """
    npts = len(x)
    # Make a vector of yerr values if yerr is a scalar.
    try:
        if len(yerr) != npts:
            raise ValueError, 'x, yerr length mismatch!'
    except TypeError:
        yerr = yerr * ones(npts, Float)
    # Check input sizes.
    if len(y) != npts:
        raise ValueError, 'x, y length mismatch!'
    if npts < 2: raise ValueError, 'Need  npts >= 2!'
    # _regress.wlss alters x, y (but not yerr) so copy them twice.
    try:
        x, y = x.copy(), y.copy()
    except AttributeError:
        x, y = array(x), array(y)
    x2, y2 = x.copy(), y.copy()
    # Handle copying yerr here so the wlss call doesn't duplicate
    # any needed conversion of a sequence to an array.
    try:
        yerr = yerr.copy()
    except AttributeError:
        yerr = array(yerr)
    inter, ierr, slope, serr, scat, error = _regress.wlss(x2, y2, yerr)
    if error == 0:
        raise ValueError, 'Input data has zero Sum(x-xbar)(y-ybar)!'
    x, y = x.copy(), y.copy()
    result = LinRegResult('wlss', slope, inter, sqrt(serr), sqrt(ierr),\
             x, y, yerr=yerr)
    # Add extra attribute.
    result.scat = sqrt(scat)
    return result


class BCES(object):
    """
    BCES:  Bivariate, Correlated Errors and Scatter
    
    Perform four linear regressions for cases where both variables
    are subject to known errors, and there is intrinsic scatter.
    Errors may be correlated or uncorrelated, and may be
    homoscedastic or heteroscedastic, i.e. the same or different
    for each datum.

    Results are available as four LinRegResult objects, stored
    as attributes:
        yonx    =   Ordinary least squares for (y|x)
        xony    =   Ordinary least squares for (x|y)
        bisect  =   Ordinary least squares bisector
        orthog  =   Orthogonal least squares
    """

    def __init__(self, x, xerr, y, yerr, covar):
        npts = len(x)
        # Make a vector of err values if scalars are provided.
        # If a vector is provided, copy it.
        try:
            if len(xerr) != npts:
                raise ValueError, 'x, xerr length mismatch!'
            try:
                xerr = xerr.copy()
            except AttributeError:
                xerr = array(xerr)
        except TypeError:
            xerr = xerr * ones(npts, Float)
        try:
            if len(yerr) != npts:
                raise ValueError, 'x, yerr length mismatch!'
            try:
                yerr = yerr.copy()
            except AttributeError:
                yerr = array(yerr)
        except TypeError:
            yerr = yerr * ones(npts, Float)
        try:
            if len(covar) != npts:
                raise ValueError, 'x, covar length mismatch!'
            try:
                covar = covar.copy()
            except AttributeError:
                covar = array(covar)
        except TypeError:
            covar = covar * ones(npts, Float)
        # Check input sizes.
        if len(y) != npts:
            raise ValueError, 'x, y length mismatch!'
        if npts < 2: raise ValueError, 'Need  npts >= 2!'
        # Handle copying x, y here so the bess call doesn't duplicate
        # any needed conversion of a sequence to an array.
        try:
            x, y = x.copy(), y.copy()
        except AttributeError:
            x, y = array(x), array(y)
        # Do the regressions.
        inter, slope, ierr, serr, serr_IFAB = \
            _regress.bess(x, xerr, y, yerr, covar)
        # Store results.
        ierr, serr, serr_IFAB = sqrt(ierr), sqrt(serr), sqrt(serr_IFAB)
        self.yonx = LinRegResult('BCES-yonx', slope[0], inter[0],
            serr[0], ierr[0], x, y, xerr, yerr)
        self.xony = LinRegResult('BCES-xony', slope[1], inter[1],
            serr[1], ierr[1], x, y, xerr, yerr)
        self.bisect = LinRegResult('BCES-bisect', slope[2], inter[2],
            serr[2], ierr[2], x, y, xerr, yerr)
        self.bisect.serr_IFAB = serr_IFAB[2]
        self.orthog = LinRegResult('BCES-orthog', slope[3], inter[3],
            serr[3], ierr[3], x, y, xerr, yerr)
        self.orthog.serr_IFAB = serr_IFAB[3]


class SixLin(object):
    """
    SixLin:  Six bivariate linear regressions
    
    Perform six linear regressions for bivariate data without
    known errors.

    Results are available as six LinRegResult objects, stored
    as attributes:
        yonx    =   Ordinary least squares for (y|x)
        xony    =   Ordinary least squares for (x|y)
        bisect  =   Ordinary least squares bisector
        orthog  =   Orthogonal least squares
        orthog  =   Orthogonal least squares
        rma     =   Reduced major axis
        mols    =   Mean ordinary least squares
    """

    def __init__(self, x, y):
        npts = len(x)
        # Check input sizes.
        if len(y) != npts:
            raise ValueError, 'x, y length mismatch!'
        if npts < 2: raise ValueError, 'Need  npts >= 2!'
        # Handle copying x, y here so the sixlin call doesn't duplicate
        # any needed conversion of a sequence to an array.
        try:
            x, y = x.copy(), y.copy()
        except AttributeError:
            x, y = array(x), array(y)
        # sixlin alters x, y so copy them twice.
        try:
            x, y = x.copy(), y.copy()
        except AttributeError:
            x, y = array(x), array(y)
        x2, y2 = x.copy(), y.copy()
        # Do the regressions.
        inter, ierr, slope, serr, error = _regress.sixlin(x2, y2)
        if error == 0:
            raise ValueError, 'Input data has zero Sum(x-xbar)(y-ybar)!'
        # Store results.
        self.yonx = LinRegResult('SixLin-yonx', slope[0], inter[0],
            serr[0], ierr[0], x, y)
        self.xony = LinRegResult('SixLin-xony', slope[1], inter[1],
            serr[1], ierr[1], x, y)
        self.bisect = LinRegResult('SixLin-bisect', slope[2], inter[2],
            serr[2], ierr[2], x, y)
        self.orthog = LinRegResult('SixLin-orthog', slope[3], inter[3],
            serr[3], ierr[3], x, y)
        self.rma = LinRegResult('SixLin-rma', slope[4], inter[4],
            serr[4], ierr[4], x, y)
        self.mols = LinRegResult('SixLin-mols', slope[5], inter[5],
            serr[5], ierr[5], x, y)
