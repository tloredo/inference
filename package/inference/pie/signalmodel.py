"""
Parameterized signal model base class.
"""

from .param import *
from .autoname import HasAutoNamed

# Exception for running out of steps:
class StopStepping(ParamError):
    """Exception for attempting to step model params beyond the last set."""
    pass

class SignalModel(HasAutoNamed):
    """
    Base class for parameterized signal models.  This must be mixed in with a
    subclass of the Inference base class that actually oversees most of the
    parameter bookkeeping.
    
    Users must override *signal* and optionally *signals*.
    """
    # ??? Should we require signal if signals is defined?

    # *** Note that none of the param state-monitoring methods can
    # be defined to use *predictor* parameters (they won't exist in the
    # user's SignalModel subclass).  Is this a
    # significant limitation?  Predictor parameters should
    # not affect calculation of the signal, so I don't think
    # it's a significant limitation.

    def signal(self, *args):
        """
        Override this method to return the value(s) defined by the model.
        
        To calculate them, access the parameter values as attributes of self
        (e.g., "self.theta").  The arguments to this method should be quantities
        *other than* the model parameters needed for evaluating the model, 
        e.g., the abscissa at which a fitted function is being evaluated.
        """
        raise NotImplementedError

    def signals(self, argvec):
        """
        Override this method to return the value(s) defined by the model for
        a vector containing many sets of arguments for the model.  Use vector
        operations, ufuncs, or custom extensions to accelerate the calculation.
        """
        raise NotImplementedError

    def on_param_change(self, param):
        """
        Note when the value of a parameter is changed.
        
        Override this if work must be done whenever a param value changes; the
        most common scenario where this may be the case is when the allowed 
        ranges of the parameters are coupled, in which case a check for
        parameter validity should be implemented via this method.  Another use
        case is when derived parameter values (that are not Param instances)
        must be updated.
        
        More common model setup chores, such as initializing interpolators,
        should be implemented via on_use, since such setup work is usually not
        necessary until actually evaluating model output, which may not happen
        until after several parameters are changed.
        
        Do *not* assign any parameters in this method; this will cause
        runaway recursive calling of `on_param_change`.
        """
        pass

    def on_use(self):
        """
        Called when the parameters of the model actually get *used* for
        calculating model signal values.  Override this, e.g., to initialize a
        table lookup or interpolator, or to do other "setup" work that need not
        change when the model value is evaluated for different arguments.
        """
        pass

    def _on_param_value_change(self, param):
        """
        A "buffer" between internal calls and the user's on_param_change
        method, insulating the user from maintaining on_use_done.
        """
        self.on_use_done = False
        self.on_param_change(param)

    def _on_use(self):
        """
        A "buffer" between internal calls and the user's on_use method,
        insulating the user from maintaining on_use_done.
        """
        self.on_use()
        self.on_use_done = True

