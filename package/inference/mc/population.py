"""
Classes and functions for describing and sampling discrete populations.
"""

# TODO:  Distinguish sizes from weights (normalized).

__author__ = "Tom Loredo"

from numpy import array, random
from _ppssampler import _ppssampler, equalprob
from _ppssampler import set_rng_state, get_rng_state
from inference.utils.pl import pl
# try:
#     import pylab as pl
# except ImportError:
#     pl = None

__all__ = ['Population', 'Population1D']

class Population(object):

    def __init__(self, items=None, weights=None):
        if weights is None:  # The equal probability case
            self.equiprob = True
            try:  # items is a sequence
                self.npopn = len(items)
                self.items = items
            except TypeError:  # items is an int
                self.npopn = int(items)
                self.items = list(range(self.npopn))
        elif items is None:  # Use indices as items
            self.equiprob = False
            self.npopn = len(weights)
            self.items = list(range(self.npopn))
            self.weights = array(weights, float)
        else:  # Use a list of items *and* weights
            self.equiprob = False
            self.npopn = len(weights)
            if len(items) != self.npopn:
                raise ValueError('Lengths of items & weights must match!')
            self.items = items
            self.weights = array(weights, float)
        
        self.did_init = False
        self.did_Sampford_init = False
        self.did_Sampford_tables = False

    def sample(self, nsamp):
        """
        Return a set of nsamp samples from the population, sampled with
        replacement.
        """
        # *** Implement equiprob case.
        if self.equiprob:
            raise NotImplementedError('Awaiting code...')
        if not self.did_init:
            self.sampler = _ppssampler(self.weights)
            self.did_init = True
        # Track the RNG state within the sampler, to update NumPy's RNG state.
        # Internally we only use the MT state; any extra state for cached
        # normal or other samples can just be copied.
        rng_state = random.get_state()
        mt_state, extra_state = rng_state[:3], rng_state[3:]
        set_rng_state(*mt_state)  # *** modify to handle full rng state
        indices = self.sampler.sample(nsamp)
        new_state = list(get_rng_state())
        new_state.extend(extra_state)
        random.set_state(new_state)
        return [self.items[i] for i in indices]

    def max_subset(self):
        """
        Return the maximum sample size for PPS sampling without replacement.
        
        The limiting size arises because PPS sampling without replacement 
        requires nsamp*(max normalized weight) <= 1.  If this is violated
        for the desired sample size, you may consider trimming the large
        weight members from the population and including them in every
        sample (of course, they will all have inclusion probability of
        unity, regardless of size).
        """
        if self.did_Sampford_init:
            return int(1./self.max_wt)
        else:
            return int(sum(self.weights)/self.weights.max())

    def subset_pps(self, nsamp):
        """
        Return a sample of nsamp distinct items from the population, sampled 
        without replacement with probability proportional to size (PPS)
        according to Sampford's sampling scheme.
        """
        # Copy the whole population if nsamp = npopn.
        if nsamp == self.npopn:
            return [item for item in self.items]
        set_rng_state(*random.get_state())
        if self.equiprob:
            pool = arange(self.npopn)
            indices = equalprob(nsamp, pool)
        else:
            # This part of setup has to be done before any sampling.
            if not self.did_init:
                print('Initing ppssampler...')
                self.sampler = _ppssampler(self.weights)
                self.did_init = True
            # This part has to be done before any sampling w/o replacement.
            if not self.did_Sampford_init:
                print('Initing wts...')
                self.sort_indices, self.sort_wts, self.tot_wt = \
                    self.sampler.prepwts(self.weights)
                self.max_wt = self.sort_wts[0]/self.tot_wt  # Max wt, normed
                self.nsamp = 0
                self.did_Sampford_init = True
                self.did_Sampford_tables = False
            # This part has to be done when sample size changes.
            if self.nsamp != nsamp:
                print('Initing ratios...')
                if nsamp > self.npopn:
                    raise ValueError('nsamp larger than population size!')
                if nsamp*self.max_wt > 1:
                    raise ValueError('Sample size too large for PPS sampling!')
                self.sampler.prepratios(nsamp, self.sort_wts, self.tot_wt)
                self.did_Sampford_tables = False
                self.nsamp = nsamp
            self.ntry, sindices = self.sampler.samplenr()
            indices = [self.sort_indices[i] for i in sindices]
        result = [self.items[i] for i in indices]
        random.set_state(get_rng_state())
        return result

    def subset_pps5(self, nsamp):
        """
        Return a sample of nsamp distinct items from the population, sampled 
        without replacement with probability proportional to size (PPS)
        according to Sampford's sampling scheme.
        
        5-table lookup samplers are used within Sampford's algorithm to
        accelerate the sampling for large populations.
        """
        # Copy the whole population if nsamp = npopn.
        if nsamp == self.npopn:
            return [item for item in self.items]
        set_rng_state(*random.get_state())
        if self.equiprob:
            pool = arange(self.npopn)
            indices = equalprob(nsamp, pool)
        else:
            # This part of setup has to be done before any sampling.
            if not self.did_init:
                print('Initing ppssampler...')
                self.sampler = _ppssampler(self.weights)
                self.did_init = True
            # This part has to be done before any sampling w/o replacement.
            if not self.did_Sampford_init:
                print('Initing wts...')
                self.sort_indices, self.sort_wts, self.tot_wt = \
                    self.sampler.prepwts(self.weights)
                self.max_wt = self.sort_wts[0]/self.tot_wt  # Max wt, normed
                self.nsamp = 0
                self.did_Sampford_init = True
            # This part has to be done when sample size changes.
            if self.nsamp != nsamp:
                print('Initing ratios...')
                if nsamp > self.npopn:
                    raise ValueError('nsamp larger than population size!')
                if nsamp*self.max_wt > 1:
                    raise ValueError('Sample size too large for PPS sampling!')
                self.sampler.prepratios(nsamp, self.sort_wts, self.tot_wt)
                self.sampler.prepratiotables()
                self.did_Sampford_tables = True
                self.nsamp = nsamp
            # This may happen if subset_pps is called before subset_pps5.
            if not self.did_Sampford_tables:
                print('Initing ratio tables...')
                self.sampler.prepratiotables()
                self.did_Sampford_tables = True
            self.ntry, indices = self.sampler.samplenr5()
            # Note the 5-table version returns unsorted indices.
            # indices = [self.sort_indices[i] for i in sindices]
        result = [self.items[i] for i in indices]
        random.set_state(get_rng_state())
        return result


class Population1D(Population):
    """
    A Population object specialized for populations indexed by a
    single (1-D) real-valued quantity that gives the "size" of
    each member.
    """

    def __init__(self, vals, weights, err=None):
        indices = array(vals).argsort()
        self.vals = vals[indices].copy()
        self.weights = array(weights)[indices].copy()
        if err == None:
            self.err = None
        else:
            self.err = array(err)[indices].copy()
        Population.__init__(self, self.vals, self.weights)
        self.cdf = self.weights.cumsum()
        self.hazard = self.cdf[::-1].copy()

    def haz_pts(self, start=None, end=None):
        """
        Return arrays of points specifying the hazard dist'n over the range
        [start, end].  The range must fully span the range of
        detected values.  Also return arrays of points specifying
        error bars, if defined on creation.
        """
        if start is None:
            start = self.vals[0]
        if end is None:
            end = self.vals[-1]
        if start>self.vals[0] or end<self.vals[-1]:
            raise ValueError('Range must span the range of sampled values!')
        # Start the descending CDF.
        absc, ord = [start], [1.]
        # Add pairs of points for each uncensored value, defining jumps.
        for x, p in zip(self.vals, self.hazard):
            absc.extend([x, x])
            ord.extend([ord[-1], p])
        # The last step is zero.
        absc.append(end)
        ord.append(0.)
        if self.err == None:
            return array(absc), array(ord)
        else:
            # For error bars, just copy the stored errors in the middle of
            # the CDF bins.
            eabsc = []
            for i in range(len(self.vals)-1):
                eabsc.append( .5*(self.vals[i]+self.vals[i+1]) )
            eabsc.append( .5*(self.vals[-1]+end) )
            return array(absc), array(ord), array(eabsc), self.hazard.copy(), \
                self.err.copy()

    def plot(self, start=None, end=None):
        """
        Plot the hazard over the range [start,end], which must span
        the range of uncensored values.
        """
        if not pl:
            raise RuntimeError('Cannot plot without pylab!')
        if start is None:
            start = self.vals[0]
        if end is None:
            end = self.vals[-1]
        if self.err == None:
            a, o = self.haz_pts(start, end)
        else:
            a, o, ea, eo, ee = self.haz_pts(start, end)
        pl.plot(a, o, 'b-', linewidth=2)
        if self.err != None:
            pl.errorbar(ea, eo, ee, fmt='o', markersize=0)
