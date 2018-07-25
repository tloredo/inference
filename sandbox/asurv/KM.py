from numpy import array, zeros, concatenate, Float
from _asurvkm import kmestm, kmdif, quart, plestm
from population import *
try:
    import pylab as pl
except ImportError:
    pl = None

__all__ = ['KM', 'plestm']


class KM(Population1D):
    """
    Kaplan-Meier estimator for right- or left-censored data.
    """

    def __init__(self, dtxns=None, ulims=None, llims=None, all=None, ind=None):
        """
        Find the Kaplan-Meier CDF estimate from right-censored (upper limits)
        or left-censored (lower limits) data.  The data may be entered as
        detected values (dtxns) and EITHER uppler (ulims) or lower (llims)
        limits, or as a single set of data values (all) accompanied by
        an array of integer indicators (0=dtxn, -1=lower limit, 1=upper limit).
        """
        if all != None:
            if dtxns or ulims or llims:
                raise ValueError, 'Illegal selection of input parameters!'
            if (ind == None) or (len(ind) != len(all)):
                raise ValueError, '"all" values must be accompanied by indicators!'
            self.all = all
            self.ind = ind
            self.ntot = len(all)
        else:
            if (ulims != None) and (llims != None):
                raise ValueError, 'Only one of upper or lower limits allowed!'
            if ulims != None:
                nlim = len(ulims)
                self.limit = 'Upper'
                self.all = concatenate([dtxns, ulims])
            elif llims != None:
                nlim = len(llims)
                self.limit = 'Lower'
                self.all = concatenate([dtxs, llims])
            else:
                raise ValueError, 'Must specify dtxns and 1 of ulims/llims!'
            nd = len(dtxns)
            self.ntot = nd + nlim
            self.ind = zeros(self.ntot)
            if self.limit == 'Upper':
                self.ind[nd:] = 1
            else:
                self.ind[nd:] = -1
        ierr, self.km_cdf, self.cdf_err, self.mean, self.std, self.nu, \
            self.uncens, self.nc, self.cens = kmestm(self.ind, self.all)
        if ierr != 0:
            raise ValueError, 'Invalid indicator array!'
        # Trim sizes of returned arrays.
        self.km_cdf = self.km_cdf[:self.nu]
        self.cdf_err = self.cdf_err[:self.nu]
        self.uncens = self.uncens[:self.nu]
        self.cens = self.cens[:self.nc]
        # Note there is no binned dist'n (yet).
        self.nbins = None
        # Define a Population1D object
        wts = zeros(len(self.uncens), Float)
        for i in range(len(self.uncens)-1):
            wts[i] = self.km_cdf[i]-self.km_cdf[i+1]
        wts[-1] = self.km_cdf[-1]
        Population1D.__init__(self, self.uncens, wts)

    def diff(self, nb, start, w):
        """
        Estimate the differential distribution in nb bins from start with
        width w.
        """
        self.nbins = nb
        self.bin_l, self.bin_u, self.bins = \
            kmdif(self.cdf, self.uncens, self.ntot, start, w, nb)

    def quartiles(self):
        """
        Return the 25%, 50%, and 75% points (quartiles) of the estimated
        distribution.  This is only valid if there are at least 4 
        uncensored values (detections).
        """
        if self.nu < 4:
            raise ValueError, 'Need > 3 uncensored values for quartiles!'
        return quart(self.uncens, self.km_cdf)

    def cdf_pts(self, start, end):
        """
        Return arrays of points specifying the KM CDF over the range
        [start, end].  The range must fully span the range of
        detected values.  Also return arrays of points specifying
        error bars.
        """
        if start>self.uncens[0] or end<self.uncens[-1]:
            raise ValueError, 'Range must span the range of uncensored values!'
        # Start the descending CDF.
        absc, ord = [start], [1.]
        # Add pairs of points for each uncensored value, defining jumps.
        for x, p in zip(self.uncens, self.km_cdf):
            absc.extend([x, x])
            ord.extend([ord[-1], p])
        # The last step is zero.
        absc.append(end)
        ord.append(0.)
        # For error bars, just copy the stored errors in the middle of
        # the CDF bins.
        eabsc = []
        for i in range(len(self.uncens)-1):
            eabsc.append( .5*(self.uncens[i]+self.uncens[i+1]) )
        eabsc.append( .5*(self.uncens[-1]+end) )
        return array(absc), array(ord), array(eabsc), self.km_cdf.copy(), \
            self.cdf_err.copy()

    def diff_pts(self, nb=None, start=None, w=None):
        """
        Return arrays of points specifying the binned distribution.
        """
        if nb is None:
            if self.nbins is None:
                raise ValueError, 'Must specify binning!'
        else:
            self.diff(nb, start, w)
        absc = 0.5*(self.bin_l + self.bin_u)
        return absc, self.bins.copy()

    def diff_histo(self, nb=None, start=None, w=None):
        """
        Return arrays of points specifying a histogram for the binned 
        distribution.
        """
        if nb is None:
            if self.nbins is None:
                raise ValueError, 'Must specify binning!'
        else:
            self.diff(nb, start, w)
        absc, ord = [self.bin_l[0]], [0.]
        for l, r, b in zip(self.bin_l, self.bin_u, self.bins):
            absc.extend([l, r])
            ord.extend([b, b])
        absc.append(absc[-1])
        ord.append(0.)
        return array(absc), array(ord)

    def plot(self, start=None, end=None):
        """
        Plot the estimated hazard funnction over the range [start,end], which 
        must span the range of uncensored values.
        """
        if not pl:
            raise RuntimeError, 'Cannot plot without pylab!'
        if start is None:
            start = min(self.vals)
        if end is None:
            end = max(self.uncens)
        a, o, ea, eo, ee = self.cdf_pts(start, end)
        pl.plot(a, o, 'b-', linewidth=2)
        pl.errorbar(ea, eo, ee, fmt='o', markersize=0)

    def plot_diff(self, nb=None, start=None, w=None):
        """
        Plot the estimated CDF over the range [start,end], which must span
        the range of uncensored values.
        """
        if not pl:
            raise RuntimeError, 'Cannot plot without pylab!'
        a, o = self.diff_histo(nb, start, w)
        pl.plot(a, o, 'g-', linewidth=2)

