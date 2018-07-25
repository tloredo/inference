from autoname import AutoNamed, HasAutoNamed
from logger import pielog

__all__ = ['Predictor', 'PredictorSet', 'DataError']


class DataError(Exception):
    pass


# IDEA: Perhaps abstract some common handler behavior in a
# PredictorHandler base class.

class PredictorHandler(object):

    def __init__(self, pinst, owner):
        """Constructor, recording the instance that owns this handler,
        and the doc string for the associated parameter.

        pinst = the Predictor instance whose attributes this handler
                will handle.
        owner = the instance of the class in which the Predictor appears.
                This will typically be an Inference subclass.
        """
        self.pinst = pinst
        self.owner = owner
        self.name = pinst.name
        self.doc = pinst.doc
        self.listeners = [owner]    # *** We don't use this capability yet.
        self._using_globals = True
        self.init()

    def init(self):
        pass

    def using_globals(self):
        return self._using_globals


class Predictor(AutoNamed):
    """An abstract base class for predictor classes that contain
    data and methods for predicting and simulating data using a signal model.

    Predictors are implemented as descriptors.  All instances of
    a predictor will share access to the same class variables storing
    the "original" values of the data.  Each instance may modify
    its working copy of the data, e.g., for simulation or sensitivity
    studies.

    Predictors may have parameters that get incorporated into the
    parameter space defining an inference (e.g., a scale factor for
    the sigmas in a Gaussian sampling distribution).

    A subclass must set self.handler_class to the handler class for
    the Predictor being defined."""

    def _report_name(self, cls, name):
        pielog.debug('Predictor instance in %s is named %s', cls, name)
        self.name = name
        self.handler_name = '__' + name + '_hndlr'
        # Keep a list of predictor names in the class dict.
        try:
            cls.predictor_names.append(name)
        except AttributeError:
            cls.predictor_names = [name]

    def clone(self):
        # ??? Can this be abstracted here?
        raise NotImplementedError

    def __get__(self, inst, owner):
        try:
            # This fails only on the 1st access.
            handler = getattr(inst, self.handler_name)
            return handler
        except AttributeError:
            # On 1st access, install the handler in the instance.
            # This will be the 1st time we know the instance to install
            # the handler in.
            handler = self.handler_class(self, inst)
            setattr(inst, self.handler_name, handler)
            return handler



class PredictorSet(HasAutoNamed):
    """A base class for classes containing one or more predictors
    together defining the sample space and sampling distributions
    for an inference."""

    def clone(self):
        """Return a PredictorSet *class* that copies the current data.
        This will be useful only if data has been simulated; this
        will allow analysis of a particular simulated data set with
        various models/methods."""
        raise NotImplementedError
