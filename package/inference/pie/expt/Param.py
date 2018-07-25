from math import exp, log
from Numeric import array, zeros, identity, Float
from NRMin import powell

################################################################################
#  Automatically named descriptors
#

"""
The following three classes comprise a type (a metaclass) and two base
classes for implementing classes with descriptors that automatically get
notified of the name of the attribute they are assigned to. The intent is
for the descriptor instance to use this name to store state in an instance
of the class containing the descriptor.  This way multiple instances
of the same descriptor can maintain different state and behavior in
a single containing class.

The motivation is to allow a class to represent a parameterized model
that may have a number of parameters of the same type (e.g., real scalars).
These parameters could not only have different values, but also different
allowed ranges (i.e., state) and different behavior (i.e., some fixed
and some varying or stepped as the model is fit to data).

Useage in a nutshell:  If a class inherting "HasAutoNamed" has descriptors
inheriting "AutoNamed", the descriptor instances will be notified of their
attribute names via calls to their "_notifyOfName" methods.
"""

class metaAutoName(type):
    """A metaclass for classes with descriptors automatically notified of thier names.

    Instances of this type can have AutoName'd descriptors that get notified of the
    names of the attributes they get assigned to."""

    def __init__(cls, name, bases, dict):
        for key,val in dict.items():
            # Search for AutoNamed descriptors.
            if isinstance(val, AutoNamed):
                # Tell the descriptor "key" is its name.
                val._notifyOfName(cls, key)
        super(metaAutoName, cls).__init__(name,bases,dict)

class AutoNamed(object):
    """An abstract base class intended for automatically passing name info
    to descriptors.  Subclasses should override "_notifyOfName" to do
    something useful with the name of the subclass instance."""

    def _notifyOfName(self, cls, name):
        """Override this method to do something useful with name."""
        raise NotImplementedError

# Note:  HasAutoNamed cannot appear before AutoNamed; the interpreter will
# complain that AutoNamed is not defined when it uses metaAutoName.

class HasAutoNamed(object):
    """Empty base class for classes with AutoNamed descriptors."""
    __metaclass__ = metaAutoName

################################################################################
#  Real-valued (i.e., float-valued) parameters
#

# Enum for parameter status:
class Status:
    undef, fixed, stepping, varying = range(4)

undef, fixed, stepping, varying = \
       Status.undef, Status.fixed, Status.stepping, Status.varying

# Enum for step types:
class StepType:
    linSteps, logSteps = range(2)

linSteps, logSteps = StepType.linSteps, StepType.logSteps

# Exceptions:
class ParamError(Exception):
    """Base class for parameter exceptions."""
    pass

class ParamRangeError(ParamError):
    """Exception for parameter out-of-range (setting or stepping)."""
    pass


class RealParamValue(float):
    """The object returned from access to a real-valued parameter.
    
    This is a Python float that forwards non-float attribute access to
    a handler that holds additional parameter state, such as the allowed
    range, the parameter nature (fixed, varying, stepped...), etc.."""

    def __new__(cls, value, handler):
        rpv = float.__new__(cls, value)
        rpv.init(handler)
        return rpv

    def __init__(self, *args, **kwds):
        pass

    def init(self, handler):
        ## print 'called RPV init', handler.name
        self.handler = handler

    def __getattr__(self, name):
        ## self.handler.printer(name)
        return getattr(self.handler, name)

class RealParamHandler(object):

    def __init__(self, name, owner, valueName, doc):
        """Constructor, recording the instance that owns this handler,
        and the doc string for the associated parameter."""
        self.name = name
        self.owner = owner  # The instance that owns this handler
        self.valueName = valueName  # Name of owner's attribute holding the value
        self.doc = doc
        self.listeners = [owner]
        self.bounded = None
        self.delta = None
        self.status = undef

    def printer(self, str):
        print 'Handler accessing attribute', str

    def addListener(self, listener):
        if not listener in self.listener:
            self.listeners.append(listener)

    def delListener(self, listener):
        try:
            self.listener.remove(listener)
        except ValueError:
            pass

    def _setStatus(self, status):
        # Note it is possible for status to change without value
        # changing, e.g. when a varied param is fixed to its current
        # value.
        old = self.status
        self.status = status
        print 'setStatus for ',self.name, old, status
        for listener in self.listeners:
            listener._onParamStatusChange(self, old, status)

    def _notifyOfStep(self):
        for listener in self.listeners:
            listener._onParamStepped(self)

    def _setValue(self, value):
        """Set the value of the parameter, without checking the range.
        This is for internal use where range checking is unnecessary."""
        self._value = RealParamValue(value, self)
        setattr(self.owner, self.valueName, self._value)
        for listener in self.listeners:
            listener._onParamValueChange(self)

    def setValue(self, value):
        """Set the value of a varying parameter."""
        if self.status != varying:
            raise ParamError, 'Parameter must be varying to set with setValue!'
        if self.bounded:
            if value>=self.lo and value<=self.hi:
                self._setValue(value)
            else:
                print 'setValue:', self.name, value
                raise ParamRangeError,"Value out of parameter's allowed range!"
        else:
            self._setValue(value)

    def getValue(self):
        if self.status == undef:
            raise ValueError, 'Parameter undefined!'
        return self._value

    def fix(self, value=None):
        """Fix the value of a parameter to the passed value, or to its
        current value if no value is passed (e.g., freeze a varying param)."""
        # *** Should this perhaps be "freeze"?
        if value == None:
            # If no value given, fix it at its current value.
            if self.status == undef:
                raise ParamError, 'Cannot fix a param with an undefined value!'
            # Note that this is a case where a param status changes
            # without the value changing.
        elif self.bounded:
            if value>=self.lo and value<=self.hi:
                self._setValue(value)
            else:
                raise ParamRangeError,"'Value out of parameter's allowed range!"
        else:
            self._setValue(value)
        self._setStatus(fixed)

    def range(self, range):
        if len(range) != 2:
            raise ValueError, 'Range must be a 2-element sequence!'
        self.lo, self.hi = float(min(range)), float(max(range))
        self.bounded = 1

    def linStep(self, lo, hi, n):
        """Step from lo to hi (inclusive) in n equal steps.
        Note: we do not check the step range vs. any existing range."""
        # *** Should I do a range check here?
        self.stype = linSteps
        self.slo, self.shi = lo, hi
        self.steps = n
        self.delta = (hi-lo)/(n-1.)
        self.stepNum = 0
        self._setValue(lo)
        self._setStatus(stepping)

    def step(self, lo, hi, n):
        """Step from lo to hi (inclusive) in n equal steps.
        Note: we do not check the step range vs. any existing range."""
        self.linStep(lo, hi, n)

    def logStep(self, lo, hi, n):
        """Step from lo to hi (inclusive) in n logarithmic steps.
        Note: we do not check the step range vs. any existing range."""
        if lo*hi <= 0.:
            raise ValueError, 'Illegal values for log stepping!'
        self.stype = logSteps
        self.slo, self.shi = lo, hi
        self.steps = n
        self.fac = exp(log(hi/lo)/(n-1.))
        self.stepNum = 0
        self._setValue(lo)
        self._setStatus(stepping)

    def nextStep(self):
        if self.status != stepping:
            raise ParamError, 'Parameter is not set to step!'
        if self.stepNum == self.steps-1:
            raise ParamRangeError, 'Requested step beyond range!'
        self.stepNum += 1
        if self.stype == linSteps:
            if self.stepNum == self.steps-1:
                self._setValue(self.shi)
            else:
                self._setValue(self._value+self.delta)
        elif self.stype == logSteps:
            if self.stepNum == self.steps-1:
                self._setValue(self.shi)
            else:
                self._setValue(self._value*self.fac)
        self._notifyOfStep()

    def prevStep(self):
        if self.status != stepping:
            raise ParamError, 'Parameter is not set to step!'
        if self.stepNum == 0:
            raise ParamRangeError, 'Requested step beyond range!'
        self.stepNum -= 1
        if self.stype == linSteps:
            if self.stepNum == 0:
                self._setValue(self.slo)
            else:
                self._setValue(self._value-self.delta)
        elif self.stype == logSteps:
            if self.stepNum == 0:
                self._setValue(self.slo)
            else:
                self._setValue(self._value/self.fac)
        self._notifyOfStep()

    def firstStep(self):
        if self.status != stepping:
            raise RuntimeError, 'Parameter is not set to step!'
        self.stepNum = 0
        self._setValue(self.slo)
        self._notifyOfStep()

    def lastStep(self):
        if self.status != stepping:
            raise RuntimeError, 'Parameter is not set to step!'
        self.stepNum = self.steps-1
        self._setValue(self.shi)
        self._notifyOfStep()

    def stepIndex(self):
        if self.status != stepping:
            raise RuntimeError, 'Parameter is not set to step!'
        return self.stepNum

    def vary(self, start=None, delta=None, range=None):
        if not range:
            if not self.bounded:
                raise ValueError, 'Need a range for varying parameter!'
        else:
            self.range(range)
        if delta != None:
            self.delta = delta
        else:
            if self.delta == None:
                raise ValueError, 'Need a delta for varying parameter!'
        self._setStatus(varying)
        # *** Can harm be done by setting status before verifying
        # the start value?
        if start != None:
            self.setValue(start)  # Note: this will verify in range.
        elif self._value == None:
            raise ValueError, 'Need a start value for varying parameter!'

    def check(self, percent):
        if self.status != varying:
            raise RuntimeError, 'Check only valid for varying parameter!'
        pc = 100.*min(self._value-self.lo, self.hi-self._value)/(self.hi-self.lo)
        return pc >= percent

    def _unboundedValue(self):
        """Return a varying param's value mapped so its range is (-inf,inf)."""
        if self.status != varying:
            raise RuntimeError, 'Unbounded access valid only for varying parameter!'
        if not self.bounded:
            raise ValueError, 'No bounds defined for parameter!'
        else:
            return log(self._value-self.lo) - log(self.hi-self._value)

    def _unboundedSet(self, uvalue):
        """Set a varying param's value via a transformed parameter that has
        its range mapped to (-inf,inf)."""
        if self.status != varying:
            raise ParamError, 'Unbounded access valid only for varying parameter!'
        # *** Is this test redundant?
        if not self.bounded:
            raise ParamError, 'No bounds defined for parameter!'
        else:
            expv = exp(uvalue)
            self.setValue( (self.lo + self.hi*expv)/(1. + expv) )
            
    def show(self):
        """Return a string containing the name, value, and doc for the param."""
        s = self.name + " = " + str(self._value)
        if self.doc:
            s += ' (' + self.doc +')'
        return s

    def addOne(self):
        """A silly test method."""
        return self._value + 1.

class RealParam(AutoNamed):
    """A real-valued (i.e., float-valued) parameter.
    
    Accessing an instance of this class returns an object that behaves like
    a float in calculations, but carries additional state and methods making
    it behave like a model parameter (e.g., ranges, parameter nature, etc.)."""

    def _notifyOfName(self, cls, name):
        print 'FloatParameter instance in',cls,'is named',name
        self.name = name
        self.valueName = '__' + name + '_val'
        self.handlerName = '__' + name + '_hndlr'
        # Keep a list of param names in the class dict.
        try:
            cls.paramNames.append(name)
        except AttributeError:
            cls.paramNames = [name]

    def __init__(self, value, doc='Undocumented float parameter.'):
        # Note: these are stored in the *class* (not inst) dict.
        self.default = value
        self.doc = doc

    def __get__(self, inst, owner):
        try:
            # This fails only if the 1st access is to the default value.
            paramValue = getattr(inst, self.valueName)
            return paramValue
        except AttributeError:
            # On 1st default access, install the default in the instance.
            handler = RealParamHandler(self.name, inst, self.valueName, self.doc)
            setattr(inst, self.handlerName, handler)
            handler.fix(self.default)
            return handler.getValue()

    def __set__(self, inst, value):
        try:
            # This will fail if we have not yet accessed the default value.
            handler = getattr(inst, self.handlerName)
        except AttributeError:
        #    raise AttributeError, 'handler should already exist!'
            handler = RealParamHandler(self.name, inst, self.valueName, self.doc)
            setattr(inst, self.handlerName, handler)
        handler.fix(value)


################################################################################
#  Parameterized model base class
#

# Exception for running out of steps:
class StopStepping(ParamError):
    """Exception for attempting to step model params beyond the last set."""
    pass

# Enum for stepping direction:
class Direction:
    next, prev = range(2)

next, prev = Direction.next, Direction.prev

def reverseDirection(drxn):
    if drxn == next:
        return prev
    elif drxn == prev:
        return next
    else:
        raise ValueError, "Invalid stepping direction!"

class ParamValues(object):
    """A container for copies of parameter values."""
    # *** If mutables are stored in a ParameterValues instance, they
    # must be explicitly copied.

    def __init__(self):
        self.paramNames = []
        self.paramDocs = {}

    def store(self, name, value, doc):
        if name in self.paramNames:
            raise ParamError, 'Parameter value already stored!'
        self.paramNames.append(name)
        self.paramDocs[name] = doc
        setattr(self, name, value)

    def __str__(self):
        s = ''
        for name in self.paramNames:
            s += name + " = " + str(getattr(self,name))
            doc = self.paramDocs[name]
            if doc:
                s += ' (' + doc +')'
            s += '\n'
        return s

class ParameterizedModel(HasAutoNamed):
    """Base class for parameterized models."""

    def __init__(self):
        self.fixedParams = []
        self.varyingParams = []
        self.steppingParams = []
        self.steppingInfo = {}  # Lists: [drxn, index, steps]
        self.onUseDone = False
        # *** Use inspect to fix the order of param names in
        # self.paramNames so show() displays params in the order
        # they were defined.

    def signal(self, *args):
        """Override this method to return the value(s) defined by the model,
        accessing the parameter values as attributes of self.  The arguments
        to this method should be quantities *other than* the model parameters
        needed for evaluating the model, e.g., the abscissa at which a fitted
        function is being evaluated."""
        raise NotImplementedError

    def _onUse(self):
        """A "buffer" between internal calls and the user's onUse method,
        insulating the user from maintaining onUseDone."""
        self.onUse
        self.onUseDone = True

    def onUse(self):
        """Called when the parameters of the model actually get *used*
        for calculating model values.  Override this, e.g., to execute
        code that derives useful constants from the parameters or to
        do other "setup" work that need not change when the model value
        is evaluated for different arguments."""
        pass

    def onParamChange(self, param):
        """Note when the value of a parameter is changed.
        Override this if work must be done whenever a param value changes; the
        most common scenario where this may be the case is when the allowed 
        ranges of the parameters are coupled, in which case a check for
        parameter validity should be implemented via this method.
        More common model setup chores should be implemented via onUse, since
        such setup work is usually required for evaluating model output,
        which may not happen until after several parameters are changed."""
        self.onUseDone = False

    def _onParamValueChange(self, param):
        """A "buffer" between internal calls and the user's onParamChange
        method, insulating the user from maintaining onUseDone."""
        self.onUseDone = False
        self.onParamChange

    def _onParamStepped(self, param):
        """Keep track of when stepping params are stepped.
        Override this if necessary."""
        # *** Eliminate this?  Is there a use case?
        pass

    def _onParamStatusChange(self, param, old, new):
        """Maintain lists of fixed, varying, and stepping variables. Initialize
        the stepping parameters whenever a parameter is changed to stepping."""
        # If old=new, don't delete and then append the param name from
        # the relevant list; this will change the order of the params
        # when the user intends only to adjust a step size, etc..
        if old == new:
            if new == stepping:
                # If a stepping param is adjusted, reset all of them.
                self.steppingInfo[param] = next
                self.resetSteps()  
            return
        if old == undef:
            if new == fixed:
                self.fixedParams.append(param)
            elif new == varying:
                self.varyingParams.append(param)
            elif new == stepping:
                self.steppingParams.append(param)
                self.steppingInfo[param] = next
                self.resetSteps()  
        if old == fixed:
            self.fixedParams.remove(param)
            if new == varying:
                self.varyingParams.append(param)
            elif new == stepping:
                self.steppingParams.append(param)
                self.steppingInfo[param] = next
                self.resetSteps()  
        if old == varying:
            self.varyingParams.remove(param)
            if new == fixed:
                self.fixedParams.append(param)
            elif new == stepping:
                self.steppingParams.append(param)
                self.steppingInfo[param] = next
                self.resetSteps()  
        if old == stepping:
            self.steppingParams.remove(param)
            self.resetSteps()
            del(self.steppingInfo[param])
            if new == fixed:
                self.fixedParams.append(param)
            elif new == varying:
                self.varyingParams.append(param)

    def resetSteps(self):
        """Reset the stepped parameters and their step directions."""
        # This is slightly wasteful: when x is set to step it is
        # initialized, and this initializes it again.  Note that
        # this will result in a duplicate call to onParamChange.
        for param in self.steppingParams:
            param.firstStep()
            self.steppingInfo[param] = next

    def nextStep(self):
        """Change the parameters to the next set on the stepped grid.
        This is done in a "zig-zag" fashion, e.g., for two parameters, (x,y),
        x is incremented up to its maximum value, then y is incremented
        *without* changing x, and subsequently x is stepped *down*.  This is
        done to facilitate calculating profiles/projections on a grid.  For
        such calculations, the current *varying* parameter values (which
        will have been optimized at the previous grid point) will only be
        a good starting point for nearby points on the grid.  If after
        incrementing y we were to reset x to its minimum value, the current
        varying parameters could be at a very bad location, greatly
        prolonging optimization or even preventing it."""
        last = self.steppingParams[-1]
        for param in self.steppingParams:
            drxn = self.steppingInfo[param]
            try:
                # Try stepping a param in its current drxn;
                # return at the 1st success.
                if drxn == next:
                    param.nextStep()
                elif drxn == prev:
                    param.prevStep()
                return
            except ParamRangeError:
                # If we bumped into the end of the last stepping param's
                # range, quit.  Otherwise, reverse the step drxn of this
                # param and try stepping the next one.
                if param == last:
                    # *** Should we resetSteps() here?  It may screw up
                    # an optimization if the user expects varying params
                    # to be in a good start state (an unlikely sit'n).
                    self.resetSteps()
                    raise StopStepping, "Attempted to step beyond the last step!"
                else:
                    self.steppingInfo[param] = reverseDirection(drxn)

    def stepNums(self):
        """Return the current indices for stepping params."""
        nums = []
        for param in self.steppingParams:
            nums.append(param.stepNum)
        return tuple(nums)

    def stepShape(self):
        """Return the numbers of steps for stepping params in a tuple.
        This would be the "shape" of an array holding results on the
        stepping grid."""
        shape = []
        for param in self.steppingParams:
            shape.append(param.steps)
        return tuple(shape)

    def _setVarying(self, *args):
        """Set the values of varying parameters.  The arguments are
        assigned to the currently varying parameters in the order they
        were declared varying."""
        if len(args) != len(self.varyingParams):
            raise ParamError, 'Wrong number of varying params!'
        for param, val in zip(self.varyingParams, args):
            param.setValue(val)

    def getParams(self):
        """Return a ParamValues instance storing the current param values
        as attributes."""
        pv = ParamValues()
        for name in self.paramNames:
            param = getattr(self ,name)
            pv.store(name, param.getValue(), param.doc)
        return pv

    def setParams(self, pv):
        """Set the values of params from a ParamValues instance.
        This is valid only if all params are either fixed or varying."""
        # *** Should this just skip stepped params rather than quit?
        # *** Should it check there are no "extras" in the passed pv?
        for name in self.paramNames:
            param = getattr(self, name)
            if param.status == fixed:
                param.fix(getattr(pv,name))
            elif param.status == varying:
                param.vary(getattr(pv,name))
            else:
                raise ParamError, 'Can only setParam when params fixed or varying!'

    def showStatus(self):
        """Print lists of fixed, stepped, and varying params."""
        f = [param.name for param in self.fixedParams]
        s = [param.name for param in self.steppingParams]
        v = [param.name for param in self.varyingParams]
        print 'Fixed params:   ', f
        print 'Stepping params:', s
        print 'Varying params: ', v
        
    def show(self):
        """Print basic parameter info: name, value, doc."""
        # *** Provide explicit support for "derived" params, or
        # expect the user to override this?  Explicit support as
        # AutoNamed descriptors may be desirable to avoid forcing
        # users to implement __init__, though derived params can
        # likely be introduced in onUse rather than __init__.
        for name in self.paramNames:
            print getattr(self, name).show()


# Exceptions for inferences:
class InferenceError(Exception):
    pass

class DataError(InferenceError):
    pass

class SampledChisqrPredictor(object):
    """Predictor for data consisting of a function sampled at known points
    and measured with noise with Gaussian uncertainties of known std dev'n."""

    # *** Should we allow parameters in Predictors as well as (signal)
    # models?  If so, we should refactor param-handling stuff in
    # ParametricModel to Inference.  We'll probably also need extra
    # code (another metaclass?) to merge the signal and predictor
    # parameter lists.

    def setData(self, *args):
        if len(args) == 3:
            self.locns, self.vals, self.sigs = args
            if len(self.locns) != len(self.vals) or \
                   len(self.vals) != len(self.sigs):
                raise DataError, 'Mismatch in lengths of arguments!'
        elif len(args) == 1:
            all = args[0]
            try:
                rows, cols = all.shape
                self.locns = all[:,0]
                self.vals = all[:,1]
                self.sigs = all[:,2]
            except AttributeError:
                # If *all* is not an array, treat it as a list.
                # Only the list format can handle >1-d locns using
                # a single setData argument.
                self.locns, self.vals, self.sigs = [], [], []
                for row in all:
                    self.locns.append(row[0])
                    self.vals.append(row[1])
                    self.sigs.append(row[2])
                self.locns = array(self.locns)
                self.vals = array(self.vals)
                self.sigs = array(self.sigs)
            except:
                raise DataError, 'Bad data array format!'
        self.ndata = len(self.vals)
        if len(self.locns.shape) == 1:
            self.dimen = 1
        elif len(self.locns.shape) ==2:
            self.dimen = self.locns.shape[1]
        else:
            raise DataError, 'Locations should be scalars or 1-d arrays!'

    def getPrdxns(self):
         # *** Add onUse calls here or in signal.
         # *** Don't duplicate signal calls if not necessary.
        try:
            self.prdxns = self.signals(self.locns)
        except:
            if self.dimen == 1:
                self.prdxns = [self.signal(locn) for locn in self.locns]
            else:
                self.prdxns = [self.signal(*locn) for locn in self.locns]
       
    def chisqr(self):
        # *** Add ability to retrieve chisqr value when already
        # calculated.
        self.getPrdxns()
        self.resids = self.vals - self.prdxns
        self._chisqr = sum((self.resids/self.sigs)**2)
        return self._chisqr

    def mockData(self, args=None, sigs=None):
        # *** Return new data, or replace original data?
        self.getPrdxns()
        # *** Sample data via Gaussians.
        raise NotImplementedError, 'Oops, unfinished!'

class SimpleValueGrid(object):
    """Simple storage for scalar results from computation on stepped param grids."""

    def __init__(self, steppingParams):
        self.ndim = len(steppingParams)
        self.info = []
        self.shape = []
        for param in steppingParams:
            self.info.append( (param.name, param.slo, param.shi, \
                               param.steps, param.stype) )
            self.shape.append(param.steps)
        self.vals = zeros(self.shape,Float)

    def __str__(self):
        """Dummy string method."""
        s = str(self.ndim) + '-d SimpleGrid instance:\n'
        for name,slo,shi,steps,stype in self.info:
            if stype == linSteps:
                s += '  %s in [%g,%g]; %d linear steps\n' % (name, slo, shi, steps)
            elif stype == logSteps:
                s += '  %s in [%g,%g]; %d log steps\n' % (name, slo, shi, steps)
        return s

class SimpleLocnGrid(object):
    """Simple storage for array results from computation on stepped param grids,
    e.g., a set of parameter estimates (of varyingParams) throughout the grid."""

    def __init__(self, steppingParams, storedParams):
        self.ndim = len(steppingParams)
        self.info = []
        self.shape = []
        for param in steppingParams:
            self.info.append( (param.name, param.slo, param.shi, \
                               param.steps, param.stype) )
            self.shape.append(param.steps)
        self.values = zeros(self.shape,Float)
        raise NotImplementedError, 'Oops, not finished yet!'

    def __str__(self):
        """Dummy string method."""
        s = str(self.ndim) + '-d SimpleGrid instance:\n'
        for name,slo,shi,steps,stype in self.info:
            if stype == linSteps:
                s += '  %s in [%g,%g]; %d linear steps\n' % (name, slo, shi, steps)
            elif stype == logSteps:
                s += '  %s in [%g,%g]; %d log steps\n' % (name, slo, shi, steps)
        return s

class Inference(object):

    def _objective(self, args):
        """The quantity maximized for inference, as a function of any
        varying parameters.  This is a kind of "curry" of the objective
        function defined by whatever predictor gets mixed in with an
        inference subclass."""
        for param, val in zip(self.varyingParams, args):
            param.setValue(val)
        ## print args, self.objective()
        return self.objective()

    def doGrid(self):
        """Evaluate the inference on the grid defined by stepped params.
        If there are varying params, the objective function will be
        maximized w.r.t. them over the grid; otherwise the objective
        will simply be evaluated."""
        if self.varyingParams:
            raise NotImplementedError, 'Optimization on grid not yet implemented!'
        results = SimpleValueGrid(self.steppingParams)
        while True:
            ind = self.stepNums()
            results.vals[ind] = self._objective([])
            try:
                self.nextStep()
            except StopStepping:
                break
        return results

    def optimize(self):
        """Optimize the objective function with respect to varying params."""
        p, d = [], []
        for param in self.varyingParams:
            p.append(param.getValue())
            d.append(param.delta)
        p = array(p)
        d = identity(len(d)) * array(d,Float)
        return powell(self._objective, p, drxns=d)

class SampledChisqrInference(Inference,SampledChisqrPredictor):

    def objective(self):
        return self.chisqr()
