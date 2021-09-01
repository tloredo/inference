from numpy import linspace, concatenate, zeros, sqrt, log, exp
from numpy import iterable, empty, any
from numpy.random import rand
from .population import Population

__all__ = ['PWLinear']

class PWLinear(object):
    """
    Use an array of (non-negative) function values on a grid to define a 
    piecewise-linear pdf.  Methods evaluate the pdf and draw samples from it.
    """

    # We'll ignore periods with marginal density below exp(smallexp)
    # times the max.  exp(-21) = 7.6e-10 (~ 1/billion).
    smallexp = -21.

    def __init__(self, lo=None, hi=None, vals=None, absc=None,
                 step='lin', logvals=False, cut=1.e-9):
        """
        Define a piecewise linear probability density from an array of 
        (non-negative) function values.  The algorithm targets use with
        very large arrays that might include large regions of negligibly
        small values and should remain efficient in that case.
        
        Parameters
        ----------
        lo, hi : real
            If provided, these define the domain of the pdf, and `vals` is
            assumed to be values on a uniform grid over the domain.  `absc`
            must not be provided if `lo` and `hi` are provided.
        absc : array_like
            If provided, this array specifies the argument (abscissa) for each 
            value in `vals`.  `lo` and `hi` must not be provided if `absc` is
            provided.  `step` is ignored.
        vals : array_like
            Non-negative values defining the pdf.  They need not be normalized.
        cut : real
            Specifies a cutoff; density values a factor of `cut` below the
            maximum density will be treated as zero.  Must be nonzero if
            `logvals` is True.
        """
        # Define the domain and abscissas.
        if not (lo is None) and not (hi is None):
            if step != 'lin':
                raise NotImplementedError('Only linear step is implemented!')
            self.lo, self.hi = lo, hi
            self.n = len(vals)
            self.step = 'lin'
            self.absc, self.delta = linspace(lo, hi, self.n, retstep=True)
        elif absc is None:
            raise ValueError('Must specify either lo and hi, or absc!')
        else:
            self.n = len(vals)
            if len(absc) != self.n:
                raise ValueError('Length mismatch between absc and vals!')
            self.lo, self.hi = absc[0], absc[-1]
            self.absc = absc
            self.step = None  # indicates arbitrary step sizes

        # Define and select the non-negligible values, 'weight values.'
        if logvals:
            wvals = vals - max(vals)
            selected = (wvals > log(cut)).nonzero()[0]
            negl = (wvals <= log(cut)).nonzero()[0]
            wvals[selected] = exp(wvals[selected])
            wvals[negl] = 0.
        else:
            wvals = vals/max(vals)
            selected = (wvals > cut).nonzero()[0]

        # Break up the selection into contiguous regions.
        # First get lists of indices to nonzero values in each region.
        regions = []
        region = [ selected[0] ]
        for n in selected[1:]:
            if n==region[-1]+1:
                region.append(n)
            else:
                regions.append(region)
                region = [n]
        regions.append(region)

        # Now adjust interior boundaries to include a bounding zero value.
        # Just keep track of the range spanned by each region.
        # Start with first region, which may start at index 0.
        lo, hi = regions[0][0], regions[0][-1]
        if lo > 0 and wvals[lo] > 0.:  # must do 2nd test in case cut=0.
            lo -= 1
        if hi < self.n-1 and wvals[hi] > 0.:
            hi += 1
        ranges = [ (lo, hi) ]
        # Run over remaining regions.
        for region in regions[1:]:
            lo, hi = region[0], region[-1]
            if wvals[lo] > 0.:
                lo -= 1
            if hi < self.n-1 and wvals[hi] > 0.:
                hi += 1
            ranges.append((lo,hi))

        # Now, by region, count intervals and calculate weights for nonzero
        # intervals.
        self.n_intvls = 0
        self.regions = []
        nz_wts = []  # wts for regions with nonzero support
        nz_intvls = []
        nz_pdf = []
        for lo,hi in ranges:
            self.n_intvls += hi-lo
            self.regions.append((lo, hi, self.absc[lo], self.absc[hi]))
            intvls = zeros((hi-lo, 2), float)
            # Keep track of the interval boundaries.
            intvls[:,0] = self.absc[lo:hi]  # lower boundaries
            intvls[:,1] = self.absc[lo+1:hi+1]  # upper boundaries
            nz_intvls.append(intvls)
            # Weights are just trapezoid areas.
            wts = 0.5*(wvals[lo:hi]+wvals[lo+1:hi+1]) * (intvls[:,1]-intvls[:,0])
            nz_wts.append(wts)
            pdf = zeros((hi-lo, 2), float)
            pdf[:,0] = wvals[lo:hi]
            pdf[:,1] = wvals[lo+1:hi+1]
            nz_pdf.append(pdf)
        self.nz_wts = concatenate(nz_wts)
        self.norm = self.nz_wts.sum()
        self.nz_wts = self.nz_wts / self.norm
        self.pdf_vals = wvals / self.norm
        self.nz_intvls = concatenate(nz_intvls)
        self.nz_pdf = concatenate(nz_pdf)
        self.nz_centers = 0.5*self.nz_intvls.sum(1)
        self.nz_widths = self.nz_intvls[:,1] - self.nz_intvls[:,0]
        self.popn = Population(weights=self.nz_wts)

    def sample(self, n=None):
        """
        Return pseudo-random samples from the piecewise linear pdf.
        
        With no argument, a single sample is returned; otherwise the
        requested number is returned as a 1-D array.
        """
        if n is None:
            index = self.popn.sample(1)[0]
            lo, hi = self.nz_intvls[index]
            plo, phi = self.nz_pdf[index]
            if plo==phi:  # handle flat intervals
                # print 'flat:', index, lo, hi, lo + (hi-lo)*rand()
                return lo + (hi-lo)*rand()
            else:
                r = (hi-lo)/(phi-plo)
                return lo - r*plo + r*sqrt(plo**2 + (phi**2 - plo**2)*rand())
                # return self.nz_centers[index]
        else:
            indices = self.popn.sample(n)
            lo, hi = self.nz_intvls[indices,0], self.nz_intvls[indices,1]
            plo, phi = self.nz_pdf[indices,0], self.nz_pdf[indices,1]
            flat = (phi==plo).nonzero()[0]  # id the flat intervals
            r = (hi-lo)/(phi-plo)  # will be NaN for flat intervals
            u = rand(n)
#             if len(flat) != 0:
#                 print plo**2 + (phi**2 - plo**2)*u
#                 print r
#                 print indices
#                 print phi-plo
            vals = lo - r*plo + r*sqrt(plo**2 + (phi**2 - plo**2)*u)
            if len(flat) != 0:
                vals[flat] = lo[flat] + (hi[flat]-lo[flat])*u[flat]
            return vals
            # return self.nz_centers[indices]

    def pdf(self, x):
        """
        Calculate the probability density for sampling at x (scalar or vector).
        
        This returns the density from the piecewise linear interpolation,
        including adjustment for the cutoff, i.e., for places with density
        below the cutoff, this returns 0.  Thus this is an accurate pdf
        for the samples returned by `sample`.
        """
        if iterable(x):
            dens = empty(len(x), float)
            for i, xval in enumerate(x):
                dens[i] = self._pdf(xval)
            return dens
        else:
            return self._pdf(x)

    def _pdf(self, x):
        """
        Calculate the probability density for sampling at a single point, x.

        This returns the density from the piecewise linear interpolation,
        including adjustment for the cutoff, i.e., for places with density
        below the cutoff, this returns 0.  Thus this is an accurate pdf
        for the samples returned by `sample`.
        """
        # Find the region containing x.
        inside = None
        for region in self.regions:
            if x >= region[2] and x <= region[3]:
                inside = region
                break
        if inside is None:  # x not in any region
            return 0.
        # Handle arbitrary step size case.
        if self.step is None:
            # Do bisection within the region to find the containing interval.
            raise NotImplementedError
        # Handle even step size case.
        elif self.step == 'lin':
            i = region[0] + int((x - region[2])/self.delta)
            return self.pdf_vals[i] + \
                (self.pdf_vals[i+1]-self.pdf_vals[i]) * \
                (x - self.absc[i])/(self.absc[i+1] - self.absc[i])
        else:  # log step case
            raise NotImplementedError('log step size not yet implemented!')
