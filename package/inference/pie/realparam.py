"""
Real-valued (i.e., float-valued) parameters.
"""

from .param import *
from scipy import exp, log, isinf, isnan

from .logger import pielog


class RealParamValue(float):
    """A float-like object returned on access to a real-valued parameter.

    This is a Python float that forwards non-float attribute access to
    a handler that holds additional parameter state, such as the allowed
    range, the parameter status (fixed, varying, stepped...), etc..

    Signature:  v = RealParamValue(value, handler)

      value = float value to be returned by "normal" value access
      handler = handler instance that will handle access to anything
                not part of a Python float
    """

    def __new__(cls, value, handler):
        rpv = float.__new__(cls, value)
        rpv.init(handler)
        return rpv

    def __init__(self, *args, **kwds):
        """As with float(), this is a no-op to preserve immutability."""
        pass

    def init(self, handler):
        self.handler = handler

    def copy_value(self):
        """
        Return the "bare" float value of this instance, for purposes
        of storing a read-only copy of the parameter value rather than the
        entire instance.
        """
        # We can't simply call copy.copy(rpv) with rpv a RealParamValue
        # instance.  copy.copy would need a handler to pass to the
        # __new__ method to create the copy.
        return float(self)

    def __getattr__(self, name):
        """Catch accesses to anything not in a normal float and pass
        them to the handler."""
        return getattr(self.handler, name)


# In the following, the names chosen for end-user methods tend to
# avoid "_" (e.g., linstep i/o lin_step), while "internal" methods more
# often use "_" for legibility.  This is to reduce typing for users
# and reflects SciPy "style" for user methods.

class RealParamHandler(ScalarParamHandler):
    """Handler that takes care of access to RealParam attributes other
    than the basic float value."""

    value_class = RealParamValue

    def init(self):
        self.bounded = None
        self.delta = None
        self.status = undef

    def set_value(self, value):
        """Set the value of a varying parameter."""
        if self.status != varying:
            raise ParamError('Parameter must be varying to set with set_value!')
        if self.bounded:
            if value >= self.lo and value <= self.hi:
                self._set_value(value)
            else:
                print('set_value -- attempted out-of-range set:', self.name, value)
                raise ParamRangeError("Value out of parameter's allowed range!")
        else:
            self._set_value(value)

    def get_value(self):
        if self.status == undef:
            raise ValueError('Parameter undefined!')
        return self._value

    def fix(self, value=None):
        """
        Fix the value of a parameter to the passed value, or to its
        current value if no value is passed (e.g., freeze a varying param).
        """
        # *** Should this perhaps be "freeze"?
        if value == None:
            # If no value given, fix it at its current value.
            if self.status == undef:
                raise ParamError('Cannot fix a param with an undefined value!')
            # Note that this is a case where a param status changes
            # without the value changing.
        elif self.bounded:
            if value >= self.lo and value <= self.hi:
                self._set_value(value)
            else:
                raise ParamRangeError("Value out of parameter's allowed range!")
        else:
            self._set_value(value)
        self._set_status(fixed)

    def range(self, *range):
        if len(range) == 2:  # Two values were passed.
            self.lo, self.hi = float(min(range)), float(max(range))
        elif len(range) == 1:  # A 2-tuple was passed.
            self.lo, self.hi = float(min(range[0])), float(max(range[0]))
        else:
            raise ValueError('Range must be a pair of scalars or a 2-element sequence!')
        self.bounded = 1

    # *** Perhaps make these so a single (unnamed) arg -> set n over current range.

    def linstep(self, lo=None, hi=None, n=10):
        """Step from lo to hi (inclusive) in n equal steps."""
        if self.bounded and lo != None:
            if lo < self.lo or lo > self.hi:
                raise ParamRangeError("lo step out of parameter's allowed range!")
        if self.bounded and hi != None:
            if hi < self.lo or hi > self.hi:
                raise ParamRangeError("hi step out of parameter's allowed range!")
        self.stype = lin_steps
        if lo is None:
            self.slo = self.lo
        else:
            self.slo = lo
        if hi is None:
            self.shi = self.hi
        else:
            self.shi = hi
        self.steps = n
        self.delta = (hi-lo)/(n-1.)
        self.step_num = 0
        self._set_value(lo)
        self._set_status(stepping)

    def step(self, lo=None, hi=None, n=10):
        """Step from lo to hi (inclusive) in n equal steps."""
        self.linstep(lo, hi, n)

    def logstep(self, lo=None, hi=None, n=10):
        """Step from lo to hi (inclusive) in n logarithmic steps."""
        if self.bounded and lo != None:
            if lo < self.lo or lo > self.hi:
                raise ParamRangeError("lo step out of parameter's allowed range!")
        if self.bounded and hi != None:
            if hi < self.lo or hi > self.hi:
                raise ParamRangeError("hi step out of parameter's allowed range!")
        if lo*hi <= 0.:
            raise ValueError('Illegal (lo,hi) for log stepping (need same sign)!')
        self.stype = log_steps
        if lo:
            self.slo = lo
        else:
            self.slo = self.lo
        if hi:
            self.shi = hi
        else:
            self.shi = self.hi
        self.steps = n
        self.fac = exp(log(hi/lo)/(n-1.))
        self.step_num = 0
        self._set_value(lo)
        self._set_status(stepping)

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
        # print('  next:', self.name, self.step_num, self._value)

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
        # print('  prev:', self.name, self.step_num, self._value)

    def first_step(self):
        if self.status != stepping:
            raise RuntimeError('Parameter is not set to step!')
        self.step_num = 0
        self._set_value(self.slo)
        # self._notify_of_step()

    def last_step(self):
        if self.status != stepping:
            raise RuntimeError('Parameter is not set to step!')
        self.step_num = self.steps-1
        self._set_value(self.shi)
        # self._notify_of_step()

    def step_index(self):
        if self.status != stepping:
            raise RuntimeError('Parameter is not set to step!')
        return self.step_num

    def vary(self, start=None, delta=None, range=None):
        if not range:
            if not self.bounded:
                raise ValueError('Need a range for varying parameter!')
        else:
            self.range(range)
        if delta is not None:
            self.delta = delta
        else:  # *** Perhaps use a default delta (& warn if so)?
            if self.delta == None:
                raise ValueError('Need a delta for varying parameter!')
        self._set_status(varying)
        # *** Can harm be done by setting status before verifying
        # the start value?
        if start != None:
            self.set_value(start)  # Note: this will verify in range.
        elif self._value == None:
            raise ValueError('Need a start value for varying parameter!')

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
            # Handle roundoff error near boundaries.
            # *** Make the constants platform-independent.
            # if self._value <= self.lo:
            # return -1.e-323
            # elif self._value >= self.hi:
            # return 1.e308
            # else:
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
                self._set_value((self.lo + self.hi*expv)/(1. + expv))

    def show(self):
        """Return a string containing the name, value, and doc for the param."""
        s = self.name + " = " + str(self._value)
        if self.doc:
            s += ' (' + self.doc + ')'
        return s


class RealParam(Param):
    """
    A real-valued (i.e., float-valued) parameter.

    Accessing an instance of this class returns an object that behaves like
    a float in calculations, but carries additional state and methods making
    it behave like a model parameter (e.g., ranges, parameter nature, etc.).

    RealParam is implemented as an autonamed descriptor (by subclassing Param).
    """

    handler_class = RealParamHandler

    def __init__(self, value, doc='Undocumented real parameter.'):
        # Note: these are stored in the *class* (not instance) dict (i.e., the
        # RealParam instance is a class variable, not an instance variable).
        # They will be copied to the instance via the handler.
        self.default = value
        self.doc = doc
