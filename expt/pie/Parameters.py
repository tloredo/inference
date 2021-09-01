from math import log, exp

# Enum for parameter status:
undef, fixed, stepping, varying = range(4)

# Enum for step types:
linSteps, logSteps = range(2)

allParams = {}

class Parameter(object):

    def __new__(cls, *args, **kwds):
        param = allParams.get(args[0])
        if param is not None:
            if len(args)==1 and len(kwds)==0:
                return param
            else:
                raise RuntimeError, 'No arguments allowed in ref to existing param!'
        param = object.__new__(cls)
        param.init(*args, **kwds)
        allParams[args[0]] = param
        return param

    def __init__(self, *args, **kwds):
        pass

    def init(self, name, default=None, desc=None):
        self.name = name
        self.dependents = []
        self.default = default
        self.__value = default  # Private to emphasize: access via methods!
        self.desc = desc
        self.bounded = None
        self.delta = None
        if self.default:
            self.setStatus(fixed)
        else:
            self.setStatus(undef)

    def addDependent(self, dependent):
        if not dependent in self.dependents:
            self.dependents.append(dependent)

    def delDependent(self, dependent):
        try:
            self.dependents.remove(dependent)
        except ValueError:
            pass

    def setStatus(self, status):
        self.status = status
        for dependent in self.dependents:
            dependent.paramStatusChange(self, status)

    def setValue(self, value):
        self.__value = value
        for dependent in self.dependents:
            dependent.paramValueChange(self, status)

    def value(self):
        if self.status == undef:
            raise ValueError, 'Parameter undefined!'
        return self.__value

    def set(self, value):
        if self.bounded:
            if value>=self.lo and value<=self.hi:
                self.__value = value
            else:
                raise ValueError,"'Value out of parameter's allowed range!"
        else:
            self.__value = value
        self.setStatus(fixed)

    def range(self, range):
        if len(range) != 2:
            raise ValueError, 'Range must be a 2-element sequence!'
        self.lo, self.hi = min(range), max(range)
        self.bounded = 1

    def linStep(self, lo, hi, n):
        """Step from lo to hi (inclusive) in n equal steps.
        Note: we do not check the step range vs. any existing range."""
        self.setStatus(stepping)
        self.stype = linSteps
        self.slo, self.shi = lo, hi
        self.steps = n
        self.delta = (hi-lo)/(n-1.)
        self.stepNum = 0
        self.__value = lo

    def logStep(self, lo, hi, n):
        """Step from lo to hi (inclusive) in n logarithmic steps.
        Note: we do not check the step range vs. any existing range."""
        if lo*hi <= 0.:
            raise ValueError, 'Illegal values for log stepping!'
        self.setStatus(stepping)
        self.stype = logSteps
        self.slo, self.shi = lo, hi
        self.steps = n
        self.fac = exp(log(hi/lo)/(n-1.))
        self.stepNum = 0
        self.__value = lo

    def nextStep(self):
        if self.status != stepping:
            raise RuntimeError, 'Parameter is not set to step!'
        self.stepNum += 1
        if self.stepNum >= self.steps:
            raise RuntimeError, 'Requested step beyond range!'
        if self.stype == linSteps:
            if self.stepNum == self.steps-1:
                self.__value = self.shi
            else:
                self.__value += self.delta
        elif self.stype == logSteps:
            if self.stepNum == self.steps-1:
                self.__value = self.shi
            else:
                self.__value *= self.fac

    def firstStep(self):
        if self.status != stepping:
            raise RuntimeError, 'Parameter is not set to step!'
        self.stepNum = 0
        self.__value = self.slo

    def vary(self, start=None, delta=None, range=None):
        if start != None:
            self.__value = start
        elif self.__value == None:
            raise ValueError, 'Need a start value for varying parameter!'
        if delta != None:
            self.delta = delta
        else:
            if self.delta == None:
                raise ValueError, 'Need a delta for varying parameter!'
        if not range:
            if not self.bounded:
                raise ValueError, 'Need a range for varying parameter!'
        else:
            self.range(range)
        self.setStatus(varying)

    def check(self, percent):
        if self.status != varying:
            raise RuntimeError, 'Check only valid for varying parameter!'
        pc = 100.*min(self.value-self.lo, self.hi-self.value)/(self.hi-self.lo)
        return pc >= percent

    def unboundedValue(self):
        if self.status != varying:
            raise RuntimeError, 'Unbounded access only valid for varying parameter!'
        if not self.bounded:
            raise ValueError, 'No bounds defined for parameter!'
        else:
            return log(self.value-self.lo) - log(self.hi-self.value)

    def unboundedSet(self, uvalue):
        if self.status != varying:
            raise RuntimeError, 'Unbounded access only valid for varying parameter!'
        if not self.bounded:
            raise ValueError, 'No bounds defined for parameter!'
        else:
            expv = exp(uvalue)
            self.__value = (self.lo + self.hi*expv) / (1. + expv)
            
    def __str__(self):
        return self.name + "=" + str(self.__value)

class Model(object):

    def __init__(self, *bareNames):
        self.params = {}
        for name in bareNames:
            self.params[name] = Parameter(name)

a = Parameter('a', 1)
b = Parameter('b')
c = Parameter('a')

print allParams
print a, b, c
