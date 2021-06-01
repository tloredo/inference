"""
param.py

Enums and exceptions for parameters and parameterized models.
Base class for parameters.
"""

from .autoname import AutoNamed
from .logger import pielog


# TODO:  Py-3.4 added an enum module in the stdlib; use it here.

# Enum for parameter status:
class Status:
    undef, fixed, stepping, varying = list(range(4))


undef, fixed, stepping, varying = \
    Status.undef, Status.fixed, Status.stepping, Status.varying

# Enum for step types:


class StepType:
    lin_steps, log_steps = list(range(2))


lin_steps, log_steps = StepType.lin_steps, StepType.log_steps


# Enum for parameter stepping direction:
class Direction:
    nextd, prevd = list(range(2))


nextd, prevd = Direction.nextd, Direction.prevd


def reverse_direction(drxn):
    """Return the enum for the reverse of the passed direction enum."""
    if drxn == nextd:
        return prevd
    elif drxn == prevd:
        return nextd
    else:
        raise ValueError("Invalid stepping direction!")


# Exceptions:
class ParamError(Exception):
    """Base class for parameter exceptions."""
    pass


class ParamRangeError(ParamError):
    """Exception for parameter out-of-range (setting or stepping)."""
    pass


class Param(AutoNamed):
    """
    A base class for parameters implemented as AutoNamed descriptors.

    This class implements the _report_name method required of AutoNamed
    descriptors.  For a parameter, this defines value and handler method names
    in the containing class allowing multiple instances to have unique state
    and behavior for the parameter.  It also maintains a list of parameter
    names in the containing class.

    This class also implements parameter attribute access.  The
    subclass must set self.handler (in __init__) to the handler class
    that will be used to handle subclass attributes.
    """

    def _report_name(self, cls, name):
        # Note that this gets called when the Param instance is
        # created, which happens when the containing *class* is
        # *defined*, not when it is instantiated.  We don't yet have
        # an *instance* of the containing class to save state in,
        # so we create and hold on to names that will hold the state
        # later, after an instance is created.
        pielog.debug('Param instance in %s is named %s', cls, name)
        self.name = name
        self.value_name = '__' + name + '_val'
        self.handler_name = '__' + name + '_hndlr'
        # Keep a list of param names in the class dict.
        try:
            cls.param_names.append(name)
        except AttributeError:
            cls.param_names = [name]

    def __get__(self, inst, owner):
        # On access (via an Inference instance), the arguments will be:
        #   self = This Param instance
        #   inst = The Inference class instance with this param
        #   owner = The Inference class with this param as a class variable
        try:
            # This fails only if the 1st access is to the default value.
            param_value = getattr(inst, self.value_name)
            return param_value
        except AttributeError:
            # On 1st default access, install the default in the instance.
            # This will be the 1st time we know the instance to install
            # the handler in.
            handler = self.handler_class(self, inst)
            setattr(inst, self.handler_name, handler)
            handler.fix(self.default)
            return handler.get_value()

    def __set__(self, inst, value):
        try:
            # This will fail if we have not yet accessed the default value.
            handler = getattr(inst, self.handler_name)
        except AttributeError:
            handler = self.handler_class(self, inst)
            setattr(inst, self.handler_name, handler)
        handler.fix(value)


# The handler base classes maintain the parameter status and value,
# and also handle broadcasts to listeners.
#
# All handlers will have listeners added and deleted the same way,
# but some types of parameters may have to maintain status and
# value differently.  E.g., for a vector param, there probably won't
# be a 'stepped' status.  So we separate out the listener list
# maintenance from status/value maintenance and broadcasting.

class ParamHandler(object):
    """
    A base class for handlers, implementing only the most
    basic common behaviors and saving basic state.

    This base class maintains the listener list.  Listener objects
    will be notified of Param events by calls to these methods:
        on_param_change, on_param_use
    Thus listeners must implement these methods.
    """

    def __init__(self, param_inst, owner):
        """
        Constructor, recording the instance that owns & accesses this handler,
        and the doc string for the associated parameter.
        *** Subclasses should not override __init__; if they do,
        they must call ParamHandler.__init__ or implement the following. ***
        """
        self.param_inst = param_inst
        self.name = param_inst.name
        self.owner = owner  # The instance that owns this handler
        self.value_name = param_inst.value_name  # Name of owner's attribute holding the value
        self.doc = param_inst.doc
        self.listeners = [owner]
        self.init()

    def init(self):
        """Subclasses should override this as necessary rather than
        overriding __init___."""
        pass

    def add_listener(self, listener):
        if not listener in self.listeners:
            self.listeners.append(listener)

    def del_listener(self, listener):
        try:
            self.listener.remove(listener)
        except ValueError:
            pass


class ScalarParamHandler(ParamHandler):
    """
    A base class for scalar-valued parameters.

    A sublcass must set self.value_class to point to the class implementing
    the returned value of the scalar parameter.
    """

    def _set_status(self, status):
        # Note it is possible for status to change without value
        # changing, e.g. when a varied param is fixed to its current
        # value.
        old = self.status
        self.status = status
        pielog.debug('setStatus for %s %s %s',self.name, old, status)
        for listener in self.listeners:
            listener._on_param_status_change(self, old, status)

# *** Is there a use case for this?
#    def _notifyOfStep(self):
#        for listener in self.listeners:
#            listener._on_param_stepped(self)

    def _set_value(self, value):
        """
        Set the value of the parameter, without checking the range.

        This is for internal use where range checking would be redundant.
        """
        self._value = self.value_class(value, self)
        setattr(self.owner, self.value_name, self._value)
        for listener in self.listeners:
            listener._on_param_value_change(self.param_inst)
