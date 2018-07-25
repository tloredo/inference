from numpy import exp, log, fromfunction
# from roots import zbracket

__all__ = ['lin_stepper', 'lin_enum_stepper', 'log_stepper', 'log_enum_stepper',
    'LogSummer', 'zbracket']

def frange(lo, hi, n):
    vals = fromfunction(lambda x: x, (n,))
    return lo+vals*(hi-lo)/(n-1.)
    
def flatten(*args):
    """Flattens the sequence of arguments so they form a single, flat sequence."""
    result = []
    for item in args:
        try:
            result += flatten(item[0])
            if item[1:]:
                result += flatten(item[1:])
        except:
            result.append(item)
    return result
    
    #print flatten(1,(2, (3,4), 5), 6)
    
def lin_stepper(lo, hi, n):
    """Step from lo to hi (inclusive) in n steps."""
    delta = (hi-lo)/(n-1.)
    for i in xrange(n):
        if i==0:
            offset = 0.
            yield lo
        elif i==n-1:
            yield hi	# Ending this way returns hi despite roundoff error.
        else:
            offset += delta
            yield lo+offset
            
def lin_enum_stepper(lo, hi, n):
    """Enumerated step from lo to hi (inclusive) in n steps."""
    delta = (hi-lo)/(n-1.)
    for i in xrange(n):
        if i==0:
            offset = 0.
            yield i, lo
        elif i==n-1:
            yield i, hi	# Ending this way returns hi despite roundoff error.
        else:
            offset += delta
            yield i, lo+offset
            
def log_stepper(lo, hi, n):
    """Step from lo to hi (inclusive) in n logarithmic steps.
	lo, hi must be > 0."""
    fac = exp(log(hi/lo)/(n-1))
    for i in xrange(n):
        if i==0:
            cur = lo
            yield cur
        elif i==n-1:
            yield hi	# End this way to return hi w/o roundoff error.
        else:
            cur *= fac
            yield cur
            
def log_enum_stepper(lo, hi, n):
    """Enumerated step from lo to hi (inclusive) in n logarithmic steps.
	lo, hi must be > 0."""
    fac = exp(log(hi/lo)/(n-1))
    for i in xrange(n):
        if i==0:
            cur = lo
            yield i, cur
        elif i==n-1:
            yield i, hi	# End this way to return hi w/o roundoff error.
        else:
            cur *= fac
            yield i, cur
            
class LogSummer:
    """
    Keep a running sum of x, with values input as log(x) and the sum
    stored as log(sum).  This prevents overflow/underflow for very large
    or very small x.  A running sum of x**2 can optionally be kept, and
    results returned as the log of the mean and std dev'n of x.
    """
    
    def __init__(self, dosig=None, start=None):
        """Initialize; can use start = (n, log(sum)) to start
		with results of a prior sum of n values.  If sig==True,
		a running sum of x**2 is also kept."""
        if start==None:
            self.n, self.logSum = 0, None
        else:
            self.n, self.logSum = start
        self.dosig = dosig
        if dosig:
            self.sigSum = None
            
    def __add__(self, val):
        self.n += 1
        if self.logSum == None:
            self.logSum = val
        elif self.logSum >= val:
            dif = val - self.logSum
            if dif > -40:	# Handles doubles w/o underflow
                self.logSum += log(1. + exp(dif))
        else:
            dif = self.logSum - val
            if dif > -40:
                self.logSum = val + log(1. + exp(dif))
            else:
                self.logSum = val
        if self.dosig:
            v2 = 2*val
            if self.sigSum == None:
                self.sigSum = v2
            elif self.sigSum >= v2:
                dif = v2 - self.sigSum
                if dif > -40:
                    self.sigSum += log(1. + exp(dif))
            else:
                dif = self.sigSum - v2
                if dif > -40:
                    self.sigSum = v2 + log(1. + exp(dif))
                else:
                    self.sigSum = v2
        return self
        
    def add_zero(self):
        """Adjust the summer for adding a zero value (log=-infty)."""
        self.n += 1
        
    def status(self):
        """Return count, log(sum), log(avg), log(sig).
		Note that sig is the std dev'n of the input values, not the
		uncertainty in the avg (which is sig/sqrt(n))."""
        if self.logSum == None:
            raise ValueError, 'Sum is zero so log is not possible!'
        ln = log(self.n)
        lmean = self.logSum - ln
        if self.dosig:
            l2 = self.sigSum - ln
            dif = 2*lmean - l2
            lsig = 0.5*l2
            if dif > -40:	# Otherwise the (-mean**2) part is below roundoff.
                lsig += 0.5*log(1. - exp(dif))
            return self.n, self.logSum, self.sigSum, lmean, lsig
        else:
            return self.n, self.logSum, lmean
            
if 0:
    test = logSummer(dosig=1)
    mn, sig = 0, 0
    l=range(1,11)
    #l.reverse()
    for n in l:
        test += log(n)
        mn += n
        sig += n**2
    print 'logSummer test:', test.status()
    mn = float(mn)/10
    from math import sqrt
    sig = sqrt(float(sig)/10 - mn**2)
    #sig = sqrt(sig/10)
    print log(mn), log(sig)
    from random import normalvariate
    test = logSummer(1)
    for n in xrange(1000):
        test += log(normalvariate(100.,3.))
    print 'logSummer with gauss:', test.status()


def zbracket(func, x1, x2, vals=None, args=None, maxtries=50, expand=1.6):
    """Bracket a zero of 1-d function func given an initial interval.
    Returns (x1, x2) where (x1, x2) is the bracketing interval.
    Raises a RuntimeError if no bracket is found.
    Optional parameters:
      vals = (func(x1), func(x2)), to save initial calls if these are
        available
      args = extra arguments for func (called with apply(func,(x,)+args))
      maxtries = maximum number of attempted expansions of the interval
      expand = expansion factor
    """
    if x1==x2:
        raise ValueError, 'zbrac requires a non-null initial interval!'
    if vals:
        f1, f2 = vals
    elif args:
        f1 = apply(func, (x1,)+args)
        f2 = apply(func, (x2,)+args)
    else:
        f1, f2 = func(x1), func(x2)
    tries = 0
    while tries < maxtries:
        tries += 1
        if f1*f2 < 0.:
            return x1, x2
        if abs(f1) < abs(f2):
            x1 += expand*(x1-x2)
            if args:
                f1 = apply(func, (x1,)+args)
            else:
                f1 = func(x1)
        else:
            x2 += expand*(x2-x1)
            if args:
                f2 = apply(func, (x2,)+args)
            else:
                f2 = func(x2)
    # We only make it here if bracketing failed.
    raise RuntimeError, 'Zero bracketing failed!'

