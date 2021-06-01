"""
Real-valued (i.e., float-valued) periodic parameters.
"""

from .param import *
from .realparam import RealParamHandler
from scipy import exp, log, isinf, isnan

from .logger import pielog

raise NotImplementedError('periodicparam.py not yet complete!')

# *** Are there situations where internally storing the modulo value will
# cause problems?  I.e., do any algorithms explicitly access a value
# and rely on differences being preserved?  Should we provide access to the
# non-modulo value?

# *** Stepping adds an increment; this could be a problem if we step across
# the period range.


class PeriodicParamHandler(RealParamHandler):
    """
    Handler that takes care of access to PeriodicParam attributes other
    than the basic float value.
    
    This is a subclass of RealParamHandler, overriding just a few methods
    (including _set_value from ScalarParamHandler) to account for parameter
    periodicity.
    """

    value_class = RealParamValue

    def init(self):
        self.bounded = None
        self.delta = None
        self.status = undef

    def _set_value(self, value):
        """
        Set the value of the parameter, without checking the range.
        This is for internal use where range checking would be redundant.
        """
        self._value = self.value_class(value, self)
        setattr(self.owner, self.value_name, self._value)
        for listener in self.listeners:
            listener._on_param_value_change(self)

    def range(self, *range):
        if len(range) == 2:  # Two values were passed.
            self.lo, self.hi = float(min(range)), float(max(range))
        elif len(range) == 1:  # A 2-tuple was passed.
            self.lo, self.hi = float(min(range[0])), float(max(range[0]))
        else:
            raise ValueError('Range must be a pair of scalars or a 2-element sequence!')
        self.bounded = 1

    def next_step(self):
        if self.status != stepping:
            raise ParamError('Parameter is not set to step!')
        if self.step_num == self.steps-1:
            raise ParamRangeError('Requested step beyond range!')
        self.step_num += 1
        if self.stype == lin_steps:
            if self.step_num == self.steps-1:
                self._set_value(self.shi)
            else:
                self._set_value(self._value+self.delta)
        elif self.stype == log_steps:
            if self.step_num == self.steps-1:
                self._set_value(self.shi)
            else:
                self._set_value(self._value*self.fac)
        # self._notify_of_step()

    def prev_step(self):
        if self.status != stepping:
            raise ParamError('Parameter is not set to step!')
        if self.step_num == 0:
            raise ParamRangeError('Requested step beyond range!')
        self.step_num -= 1
        if self.stype == lin_steps:
            if self.step_num == 0:
                self._set_value(self.slo)
            else:
                self._set_value(self._value-self.delta)
        elif self.stype == log_steps:
            if self.step_num == 0:
                self._set_value(self.slo)
            else:
                self._set_value(self._value/self.fac)
        # self._notify_of_step()

    def check(self, percent=1.):
        if self.status != varying:
            raise RuntimeError('Check only valid for varying parameter!')
        pc = 100.*min(self._value-self.lo, self.hi-self._value)/(self.hi-self.lo)
        return pc >= percent

    def _unbounded_get_value(self):
        """Return a varying param's value mapped so its range is (-inf,inf)."""
        if self.status == undef:
            raise ValueError('Parameter undefined!')
        if self.status != varying:
            raise RuntimeError('Unbounded access valid only for varying parameter!')
        if not self.bounded:
            raise ValueError('No bounds defined for parameter!')
        else:
##             # Handle roundoff error near boundaries.
##             # *** Make the constants platform-independent.
##             if self._value <= self.lo:
##                 return -1.e-323
##             elif self._value >= self.hi:
##                 return 1.e308
##             else:
            return log(self._value-self.lo) - log(self.hi-self._value)

    def _unbounded_set_value(self, uvalue):
        """Set a varying param's value via a transformed parameter that has
        its range mapped to (-inf,inf)."""
        if self.status != varying:
            raise ParamError('Unbounded access valid only for varying parameter!')
        # *** Is this test redundant?
        if not self.bounded:
            raise ParamError('No bounds defined for parameter!')
        else:
            expv = exp(uvalue)
            if isinf(expv) or isnan(expv):
                if uvalue > 0:
                    self._set_value(self.hi)
                else:
                    self._set_value(self.lo)
            else:
                # Use _set_value, since this should be in range.  In fact,
                # I've found just-barely out-of-range results (presumably due
                # to roundoff at limits) that trigger exceptions using set_value.
                self._set_value( (self.lo + self.hi*expv)/(1. + expv) )
            
    def show(self):
        """Return a string containing the name, value, and doc for the param."""
        s = self.name + " = " + str(self._value)
        if self.doc:
            s += ' (' + self.doc +')'
        return s


class PeriodicParam(Param):
    """
    A real-valued (i.e., float-valued) periodic parameter.
    
    Accessing an instance of this class returns an object that behaves like
    a float in calculations, but carries additional state and methods making
    it behave like a model parameter (e.g., ranges, parameter nature, etc.).

    This is implemented as an autonamed descriptor (by subclassing Param).
    """

    handler_class = PeriodicParamHandler

    def __init__(self, value, lo, hi, doc='Undocumented periodic parameter.'):
        # Note: these are stored in the *class* (not instance) dict (i.e., the
        # RealParam instance is a class variable, not an instance variable).
        # They will be copied to the instance via the handler.
        self.default = value
        self.intvl_lo = lo
        self.intvl_hi = hi
        self.period = hi - lo
        self.doc = doc
