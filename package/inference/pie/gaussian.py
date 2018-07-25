from predictor import Predictor, PredictorHandler, DataError
from numpy import array, concatenate
import scipy.stats

# A standard normal dist'n (0 mean, unit variance).
# stdnorm.rvs() returns a sample; stdnorm.rvs(n) returns n samples.
stdnorm = scipy.stats.norm(loc=0, scale=1)


# *** Perhaps this should be "GaussianHandler" with a virtual get_prdxns &
# set_data methods.  Binned and channel Predictors should be able to reuse
# everything else.  Or must they all subclass PredictorHandler directly?

class SampledGaussianHandler(PredictorHandler):
    """The handler that handles the method lookups for instances of a Predictor.

    An instance of this class is returned when the Predictor is accessed.
    Thus methods and attributes of this class appear as if they were
    methods and attributes of the Predictor.  This implements instance-specific
    state and behavior, even though the Predictor is a class variable."""
    
    def init(self):
        """
        Initialize data to point to the global data.
        """
        # Make local references to the 'global' data.
        # If we ever modify the data (for simulation or sensitivity
        # analysis), we'll make sure to not touch the original data.
        # For 'templates' we should just copy things right away.
        self.locns, self.vals = self.pinst._locns, self.pinst._vals
        self.sigmas = self.pinst._sigmas
        self.ndata, self.dimen = self.pinst.ndata, self.pinst.dimen
        self._using_globals = True  # ??? Do this here or in __init__?
        
    def localize_data(self):
        """
        Change the instance's data to be copies of the current data rather
	    than references to it (e.g., to the original data stored in class
	    rather than instance variables).
	    """
        # Note this presumes the data have a copy method (as do numpy
        # arrays).  For other data types we might need to use the copy
        # or deepcopy modules.
        if self._using_globals:
            self.locns = self.locns.copy()
            self.vals = self.vals.copy()
            self.sigmas = self.sigmas.copy()
            self._using_globals = False
        
    def changing_data(self):
        """
        Note an operation that is changing the values of the data.
        This means we need to localize the data (to avoid affecting
        other instances of this Predictor), and also make sure we
        clear any stored predictions.
        """
        # *** Keep track of each prediction separately, so any
        # re-evaluations calculate predictions only for changed data.
        self.localize_data()
        # self._data_changed  # ??? Do we need a global?
        self.prdxns = None
        
    def get_prdxns(self):
    # >>>>>> This is the ONLY place where signal() or signals()
    # >>>>>> should be explicitly called!
    # If calls are made anywhere else, make sure _on_use() gets
    # called beforehand.
    # *** Don't duplicate signal calls if not necessary.
        self.owner._on_use()  # ??? Should this be in SignalModel (monitoring param changes)?
        try:
            self.prdxns = self.owner.signals(self.locns)
        except NotImplementedError:
            if self.dimen == 1:
                self.prdxns = [self.owner.signal(locn) for locn in self.locns]
            else:
                self.prdxns = [self.owner.signal(*locn) for locn in self.locns]
                
    def log_like(self):
    # *** Memoize:  Add ability to retrieve log_like value when already
    # calculated.
        self.get_prdxns()
        self.resids = self.vals - self.prdxns
        self._log_like = -0.5*sum((self.resids/self.sigmas)**2)
        return self._log_like
        
    def chisqr(self):
    # *** Memoize:  Add ability to retrieve chisqr value when already
    # calculated.
        self.get_prdxns()
        self.resids = self.vals - self.prdxns
        self._chisqr = sum((self.resids/self.sigmas)**2)
        return self._chisqr
        
    def sim(self, locns=None, sigmas=None):  # ??? Switch based on args?
        """
        Replace existing data with simulated values.
        
        With no arguments, new data values are simulated using the existing
        locations and sigmas; they will replace the previous values.
        
        Note that the new data will be local to self (i.e., they will not
        overwrite the global values defined in the predictor).
        
        If either locns or sigmas is present, the stored locns and/or
        sigmas will be replaced with the passed values before simulation.
        
        If only one of them is present, its length must match the current
        number of data.  If both are present, the number of data can be
        different than that in the original data set.
        """
        self.changing_data()
        # Change the data locns/sigmas if necessary.
        if locns and sigmas:
            if len(locns) != len(sigmas):
                raise ValueError, 'locns/sigmas length mismatch!'
            self.locns = locns
            self.sigmas = sigmas
        elif locns:
            if len(locns) != self.ndata:
                raise ValueError, 'Incorrect number of locns!'
            self.locns = locns
        elif sigmas:
            if len(sigmas) != self.ndata:
                raise ValueError, 'Incorrect number of sigmas!'
            self.sigmas = sigmas
        self.get_prdxns()
        # Sample data via Gaussians.
        self.vals = self.prdxns + self.sigmas*stdnorm.rvs(self.ndata)
    
    # *** sim1, add_sim only simulate the NEW data.  But params could
    # have changed since the last simulation, in which case there will be an
    # internally inconsistency among the data.  We should test for this and
    # raise an exception.
    
    def sim1(self, n, locn=None, sigma=None):
        """
        Replace a single datum with a simulated value, optionally altering the
        location and/or sigma for the datum.
        """
        self.changing_data()
        if locn:
            # *** Verify locn dimension!
            self.locns[n] = locn
        if sigma:
            self.sigmas[n] = sigma
        self.get_prdxns()
        # Sample datum via Gaussian.
        self.vals[n] = self.prdxns[n] + self.sigmas[n]*stdnorm.rvs()
        
    def add_sim(self, locn=None, sigma=None, locns=None, sigmas=None):
        """
        Simulate *additional* data.
        
        If two arguments are passed (locn, sigma), a single new datum is
        simulated with the given location, sigma.
        
        If arrays are passed as locns & sigmas, simulate the appropriate
        mutiple additional data.
        """
        self.changing_data()
        # First adjust locns and sigmas.
        if locn and sigma:
            self.ndata += 1
            self.locns = concatenate( (self.locns, (locn,)) )
            self.sigmas = concatenate( (self.sigmas, (sigma,)) )
        elif locns and sigmas:
            if len(locns) != len(sigmas):
                raise ValueError, 'locns/sigmas length mismatch!'
            self.ndata += len(locns)
            self.locns = concatenate( (self.locns, locns) )
            self.sigmas = concatenate( (self.sigmas, sigmas) )
        self.ndata += len(self.locns)
        if len(self.locns.shape) == 1:
            self.dimen = 1
        elif len(self._locns.shape) == 2:
            self.dimen = self.locns.shape[1]
        else:
            raise DataError, 'Locations should be scalars or 1-d arrays!'
        # Update predictions.
        self.get_prdxns()
        

# *** Add method to restore to original (not simulated) data?
# May have to allow explicit caching to define "original" data
# for simulation templates.


    ## These are failed attempts to prevent user access to global data.
    ## Probably should make the data hidden.
    def __getitem__(self, key):
        print 'Attempting item access of', key
        
    def __setitem__(self, key, value):
        print 'Attempting item change of', key
        
        ##    def __getattribute__(self, name):
        ##        print 'Attempting direct attribute access for', name
        ##        if name in ['locns', 'vals', 'sigmas'] \
        ##                and object.__getattribute__(self, 'usingGlobals'):
        ##           raise RuntimeError, 'Use .localize_data() before manipulating data!'
        ##        return object.__getattribute__(self, name)
        
        
class SampledGaussian(Predictor):
    """
    Predictor for data consisting of a function sampled at fixed, known points
    and measured with noise with Gaussian uncertainties of known std dev'n.
    """
    
    # This is implemented as a descriptor; its instances are meant to appear as
    # *class* variables in Inference or PredictorSet (sub)classes.  This means
    # that attributes of 'self' are stored in those classes' dictionaries, not
    # in the dict of an instance of such classes.  Thus we use self only for
    # storing copies of the 'original' data associated with the predictor;
    # this becomes 'global' (i.e., cross-instance) data accessible to every
    # instance of a given Predictor. Instance-specific state (modified data,
    # predictions, etc.) is maintained in the instance dict via a handler.
    #
    # To emphasize that the original data is accessible to multiple
    # instances and thus should be considered 'immutable,' such data are
    # stored in attributes with leading underscores.
    
    handler_class = SampledGaussianHandler
    
    def __init__(self, *args, **kwds):
        print 'args:', args
        try:
            self.doc = kwds['doc']
        except KeyError:
            self.doc = None
        if args:
            self.set_data(*args)
        else:
        # *** Implement this; need locns, sigmas; or require explicit values=None.
        # This will be for pure simulation studies.
        # !!! This string confuses Eclipse's PyDev syntax colorizer for later strings.
            raise NotImplementedError('"Template" predictors not yet implemented!')
            
    def set_data(self, *args):
        """
        Set the data from arrays of locns, values & sigmas or an array of triples.
        
        If three arguments are passed, they are interpreted as arrays of
        locations, values, and sigmas respectively.
        
        If a single array is passed, it is interpreted as an array of 3-tuples,
        with each tuple of the form (locn, value, sigma) for a datum.
        """
        if len(args) == 3:
            self._locns, self._vals, self._sigmas = args
            if len(self._locns) != len(self._vals) or \
            len(self._vals) != len(self._sigmas):
                raise DataError, 'Mismatch in lengths of arguments!'
        elif len(args) == 1:
            all = args[0]
            try:
                rows, cols = all.shape
                self.locns = all[:,0]
                self.vals = all[:,1]
                self.sigmas = all[:,2]
            except AttributeError:
            # If *all* is not an array, treat it as a list.
            # Only the list format can handle >1-d locns using
            # a single set_data argument.
            # *** Not really - the dimension is shape[1]-2.
                self._locns, self._vals, self._sigmas = [], [], []
                for row in all:
                    self._locns.append(row[0])
                    self._vals.append(row[1])
                    self._sigmas.append(row[2])
                self._locns = array(self.locns)
                self._vals = array(self.vals)
                self._sigmas = array(self.sigmas)
            except:
                raise DataError, 'Bad data array format!'
                # *** Should ndata, dimen be _ndata...?  Or eliminate _...?
                # Lean toward latter, since handler will catch accesses (test this!).
        self.ndata = len(self._vals)
        if len(self._locns.shape) == 1:
            self.dimen = 1
        elif len(self._locns.shape) == 2:
            self.dimen = self.locns.shape[1]
        else:
            raise DataError, 'Locations should be scalars or 1-d arrays!'
