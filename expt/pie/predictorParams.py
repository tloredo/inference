from pie import *
from pie.predictor import *
from pie.param import Param

class TestPredictor(Predictor):

    def __init__(self, doc=None, sigfac=1., *args, **dict):
        print 'args:', args
        self.doc = doc
        self.sigfac = sigfac

    def __get__(self, inst, owner):
        try:
            # This fails only on the 1st access.
            handler = getattr(inst, self.handler_name)
            return handler
        except AttributeError:
            # On 1st default access, install the default in the instance.
            # This will be the 1st time we know the instance to install
            # the handler in.
            handler = TestHandler(self, inst)
            setattr(inst, self.handler_name, handler)
            return handler

class TestHandler0(object):

    def __init__(self, pinst, owner):
        """Constructor, recording the instance that owns this handler,
        and the doc string for the associated parameter.

        pinst = the Predictor instance whose attributes this handler
                will handle.
        owner = the instance of the class in which the Predictor appears."""
        self.pinst = pinst
        self.owner = owner
        self.name = pinst.name
        self.doc = pinst.doc
        if pinst.sigfac is None:
            self.sigfac = 1.
        else:
            getattr(owner,pinst.sigfac.name).add_listener(self)
            self.sigfac = getattr(owner, pinst.sigfac.name)
        self.listeners = [owner]    # *** We don't use this capability yet.
        # Make local references to the 'global' data.
        # If we ever modify the data (for simulation or sensitivity
        # analysis), we'll make sure to not touch the original data.
        # For 'templates' we should just copy things right away.

class TestHandler(object):

    def __init__(self, pinst, owner):
        """Constructor, recording the instance that owns this handler,
        and the doc string for the associated parameter.

        pinst = the Predictor instance whose attributes this handler
                will handle.
        owner = the instance of the class in which the Predictor appears."""
        self.pinst = pinst
        self.owner = owner
        self.name = pinst.name
        self.doc = pinst.doc
        self.listeners = [owner]    # *** We don't use this capability yet.
        self._set_predictor_params()
        # Make local references to the 'global' data.
        # If we ever modify the data (for simulation or sensitivity
        # analysis), we'll make sure to not touch the original data.
        # For 'templates' we should just copy things right away.

    def _set_predictor_params(self):
        """Called by the constructor, this finds references to parameters in the
        Predictor instance served by this handler, collects the reference names,
        and makes this handler a listener for each parameter, catching value
        changes."""
        # *** Is there any trouble if (unlikely!) a Predictor has
        # multiple refs to the same param?  Is there a use case?
        self. ppnames = []
        for key, val in self.pinst.__dict__.items():
            if isinstance(val, Param):
                self.ppnames.append(key)
                param_name = self.pinst.__dict__[key].name
                getattr(self.owner, param_name).add_listener(self)
                self.__dict__[key] = getattr(self.owner, param_name)
        pielog.debug('Params passed to this Predictor: %s',self.ppnames)

    def _on_param_status_change(self, param, old, new):
        """Called by any predictor parameters we're listening to.
        A no-op since we don't care about status changes, only value changes."""
        pass

    def _on_param_value_change(self, param):
        """Called by any predictor parameters we're listening to.
        This updates the local copies of the parameters."""
        # *** This version is slightly wasteful; it updates
        # every param, rather than just the one that was changed.
        for key in self.ppnames:
            # Find the name of the actual param instance this
            # predictor param is referring to.
            param_name = self.pinst.__dict__[key].name
            # Copy the param's value.
            self.__dict__[key] = getattr(self.owner, param_name)

    def __getattr__(self, name):
        """Catch accesses to predictor parameters left at default values
        by passing the access up to the predictor instance dict."""
        return self.pinst.__dict__[name]

class MyModel(SignalModel):

    p = RealParam(3.4)

    def signal(self, x):
        return 0.

class MyData(PredictorSet):

    q = RealParam(1.2)
    d1 = TestPredictor(doc='d1', sigfac=q)

    d2 = TestPredictor(doc='d2')

    r = RealParam(5.6)
    d3 = TestPredictor(doc='d3', sigfac=r)

class MyInf(ChisqrInference, MyModel, MyData):
    pass

inf = MyInf()

# if _init_ has self.sigfac = pinst.sigfac,
# then inf.d1.sigfac = MyData.__dict__['q']
#     = <pie.realparam.RealParam object at 0x4fefb0>

# inf.d1.sigfac.__get__(inf,type(inf)) returns the correct value,
# even with multiple MyInf instances.
