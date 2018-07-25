"""
Utilities and base class for inferences.
"""

from sets import Set
import inspect

from param import *
from autoname import HasAutoNamed
from signalmodel import SignalModel, StopStepping
from scipy import zeros, array, identity
from NRMin import powell
from predictor import PredictorSet
from containers import ParamValues, ScalarGrid, VectorGrid

# TODO:  Work out the minimizer situation; powell is imported from my
# NRMin module (not in inference), based on Numerical Recipes.  Rewrite
# it for inference, or use something in scipy.optimize.

# Enum for extremize direction:
class ExtremizeDrxn:
    maximize, minimize = 'max', 'min'

maximize, minimize = ExtremizeDrxn.maximize, ExtremizeDrxn.minimize

# Exceptions for inferences:

class InferenceError(Exception):
    pass


class Inference(HasAutoNamed):
    """
    Base class for inferences, managing basic parameter bookkeeping:
    fixing, stepping, and varying parameters.

    A subclass of Inference will typically get mixed in with a SignalModel
    subclass and a PredictorSet subclass, or a ProbModel subclass.  Inference
    thus relies on methods present in these mixin classes.

    Note:  If the inference inheriting this class defines an __init__ method,
    it should call Inference.__init__() within that method.
    """

    def __init__(self):
        self.fixed_params = []
        self.varying_params = []
        self.stepping_params = []
        self.stepping_info = {}  # vals = lists: [drxn, index, steps]
        # This sets whether to use unbounded transforms of varying params.
        self.use_unbounded = False
        self.on_use_done = False

        # ??? Keep lists of current values of various param types?

        # Default minimization method & tolerance: if a subclass
        # hasn't defined these, set them to None.
        # We don't provide a default method & tolerance here
        # because different inference methods may have different
        # appropriate default methods.
        if not hasattr(self, 'min_method'):
            self.min_method = None
        if not hasattr(self, 'min_tol'):
            self.min_tol = None

        # log_dens indicates whether the objective function can
        # be interpreted as the negative log of a distribution
        # on parameter space.  This gets used to determine if
        # results on grids can be meaningfully integrated.
        #if not hasattr(self, 'minus_log_df'):
        #    self.log_dens = False
        if not hasattr(self, 'minus_log_df'):
            self.log_dens = False

        # did_signals keeps track of whether we've already evaluated
        # the model.
        # *** NOT YET IMPLEMENTED
        # self.did_signals = False

        # *** Provide hooks for storing state or calculating candidate state
        # before storing it, to facilitate MCMC.

        # Reorder the param names associated with this inference to reflect
        # the order in the source code (the metaclass will order them in a
        # shuffled order, the order of keys in the class dict).
        # First, get the SignalModel and Predictors mixin classes.
        for cls in self.__class__.__bases__:
            if SignalModel in cls.__bases__:
                self.signal_class = cls
            elif PredictorSet in cls.__bases__:
                self.pred_class = cls
        # Get the names of params and predictors from the mixins.
        # Get them from source if available, to preserve the order
        # in which the user defined them.
        # Be careful gathering params; the PredictorSet mixin class
        # may have params.
        try:
            unordered = self.signal_class.param_names[:]
            self.params = self._get_signal_names(unordered)
        except:
            self.params = []
        try:
            unordered = self.pred_class.param_names[:]
        except:
            unordered = []
        oparams, self.predictors = self._get_pred_names(unordered)
        self.params.extend(oparams)
        # Make sure params are unique; i.e., that the Predictors class
        # hasn't duplicated a param name.
        if len(Set(self.params)) != len(self.params):
            raise InferenceError, 'A parameter name is duplicated within this inference!'

    def _get_signal_names(self, signal_params):
        """Use the source for the SignalModel and Predictors mixins to
        fix the order of param and predictor names.  If source is
        unavailable, just grab the names from the mixin class dicts."""
        # Be careful gathering params; the Predictors mixin class
        # may have params.
        oparams = []
        try:
            source = inspect.getsourcelines(self.signal_class)[0]
            for line in source:
                line = line.strip()
                for param in signal_params:
                    # *** Use a regexp here to match any white space b/4 '='.
                    if line.startswith(param+' =') or line.startswith(param+'='):
                        oparams.append(param)
                        signal_params.remove(param)
            if signal_params != []:
                raise InferenceError, 'Some params not defined in SignalModel source!'
        except IOError:  # Only in the unlikely event of interactive model def'n
            print 'Could not locate SignalModel source; proceeding with unordered params.'
            oparams = self.signal_class.param_names
        except TypeError:  # This is needed due to an apparent bug in inspect (failure can raise TypeError).
            print 'Could not locate SignalModel source; proceeding with unordered params.'
            print '(Note: This may be due to a bug in inspect.py)'
            oparams = self.signal_class.param_names
        return oparams

    def _get_pred_names(self, pred_params):
        """Use the source for the SignalModel and Predictors mixins to
        fix the order of param and predictor names.  If source is
        unavailable, just grab the names from the mixin class dicts."""
        # Be careful gathering params; the Predictors mixin class
        # may have params.
        oparams = []
        try:
            source = inspect.getsourcelines(self.pred_class)[0]
            if pred_params:
                for line in source:
                    line = line.strip()
                    for param in pred_params:
                        if line.startswith(param+' =') or line.startswith(param+'='):
                            oparams.append(param)
                            pred_params.remove(param)
                if pred_params != []:
                    raise InferenceError, 'Some params not defined in Predictors source!'
            # Now do the same to get the Predictor order.
            opreds = []
            preds = self.pred_class.predictor_names[:]
            for line in source:
                line = line.strip()
                for pred in preds:
                    if line.startswith(pred+' =') or line.startswith(pred+'='):
                        opreds.append(pred)
                        preds.remove(pred)
            if preds != []:
                raise InferenceError, 'Some Predictor(s) not defined in Predictors source!'
        except IOError:  # Only in the unlikely event of interactive pred. def'n
            print 'Could not locate PredictorSet source; proceeding with unordered params.'
            oparams.extend(self.pred_class.param_names)
            opreds = self.pred_class.predictor_names
        except:  # This is needed due to an apparent bug in inspect (failure can raise TypeError).
            print 'Could not locate PredictorSet source; proceeding with unordered params.'
            oparams.extend(self.pred_class.param_names)
            opreds = self.pred_class.predictor_names
        return oparams, opreds
        # Make sure params are unique; i.e., that the Predictors class
        # hasn't duplicated a param name.
        if len(Set(self.params)) != len(self.params):
            raise InferenceError, 'A parameter name is duplicated within this inference!'

    def _on_param_status_change(self, param, old, new):
        """Maintain lists of fixed, varying, and stepping variables. Initialize
        the stepping parameters whenever a parameter is changed to stepping."""
        # If old=new, don't delete and then append the param name from
        # the relevant list; this will change the order of the params
        # when the user intends only to adjust a step size, etc..
        if old == new:
            if new == stepping:
                # If a stepping param is adjusted, reset all of them.
                self.stepping_info[param] = next
                self.reset_steps()
            return
        if old == undef:
            if new == fixed:
                self.fixed_params.append(param)
            elif new == varying:
                self.varying_params.append(param)
            elif new == stepping:
                self.stepping_params.append(param)
                self.stepping_info[param] = next
                self.reset_steps()
        if old == fixed:
            self.fixed_params.remove(param)
            if new == varying:
                self.varying_params.append(param)
            elif new == stepping:
                self.stepping_params.append(param)
                self.stepping_info[param] = next
                self.reset_steps()
        if old == varying:
            self.varying_params.remove(param)
            if new == fixed:
                self.fixed_params.append(param)
            elif new == stepping:
                self.stepping_params.append(param)
                self.stepping_info[param] = next
                self.reset_steps()
        if old == stepping:
            self.stepping_params.remove(param)
            self.reset_steps()
            del(self.stepping_info[param])
            if new == fixed:
                self.fixed_params.append(param)
            elif new == varying:
                self.varying_params.append(param)

    def map_bounds(self):
        """
        Use transformed versions of varying parameters in optimizations, mapping
        the parameter bounds to infinity so the mapped parameter space is
        unbounded.

        This is an unsophisticated way to enable optimization with basic
        optimizers when parameters have fixed parameter boundaries.
        """
        self.use_unbounded = True

    def unmap_bounds(self):
        """
        Switch back to using the "true" (bounded) values of varying parameters
        in optimization.
        """
        self.use_unbounded = False

    def reset_steps(self):
        """Reset the stepped parameters and their step directions."""
        # This is slightly wasteful: when x is set to step it is
        # initialized, and this initializes it again.  Note that
        # this will result in a duplicate call to on_param_change.
        for param in self.stepping_params:
            param.first_step()
            self.stepping_info[param] = next

    def next_step(self):
        """
        Change the parameters to the next set on the stepped grid.

        This is done in a "zig-zag" fashion, e.g., for two parameters, (x,y),
        x is incremented up to its maximum value, then y is incremented once
        *without* changing x, and subsequently x is stepped *down*.  This is
        done to facilitate calculating profiles/projections on a grid.  For
        such calculations, the current *varying* parameter values (which
        will have been optimized at the previous grid point) will only be
        a good starting point for nearby points on the grid.  If after
        incrementing y we were to reset x to its minimum value, the current
        varying parameters could be at a very bad location, greatly
        prolonging optimization or even preventing it.
        """
        last = self.stepping_params[-1]
        for param in self.stepping_params:
            drxn = self.stepping_info[param]
            try:
                # Try stepping a param in its current drxn;
                # return at the 1st success.
                if drxn == next:
                    param.next_step()
                elif drxn == prev:
                    param.prev_step()
                return
            except ParamRangeError:
                # If we bumped into the end of the last stepping param's
                # range, quit.  Otherwise, reverse the step drxn of this
                # param and try stepping the next one.
                if param == last:
                    # *** Should we reset_steps() here?  It may screw up
                    # an optimization if the user expects varying params
                    # to be in a good start state (an unlikely sit'n).
                    self.reset_steps()
                    raise StopStepping, "Attempted to step beyond the last step!"
                else:
                    self.stepping_info[param] = reverse_direction(drxn)

    def step_nums(self):
        """Return the current indices for stepping params."""
        nums = []
        for param in self.stepping_params:
            nums.append(param.step_num)
        return tuple(nums)

    def step_shape(self):
        """Return the numbers of steps for stepping params in a tuple.
        This would be the "shape" of an array holding results on the
        stepping grid."""
        shape = []
        for param in self.stepping_params:
            shape.append(param.steps)
        return tuple(shape)

    def _set_varying(self, *args):
        """Set the values of varying parameters.  The arguments are
        assigned to the currently varying parameters in the order they
        were declared varying."""
        if len(args) != len(self.varying_params):
            raise ParamError, 'Wrong number of varying params!'
        for param, val in zip(self.varying_params, args):
            param.set_value(val)

    def get_params(self):
        """Return a ParamValues instance storing the current param values
        as attributes."""
        pv = ParamValues()
        for name in self.params:
            param = getattr(self, name)
            # pv.store(name, param.get_value(), param.doc)
            pv.store(name, param)
        return pv

    def set_params(self, pv):
        """Set the values of params from a ParamValues instance.
        This is valid only if all params are either fixed or varying."""
        # *** Should this just skip stepped params rather than quit?
        # *** Should it check there are no "extras" in the passed pv?
        for name in self.params:
            param = getattr(self, name)
            if param.status == fixed:
                param.fix(getattr(pv,name))
            elif param.status == varying:
                param.vary(getattr(pv,name))
            else:
                raise ParamError, 'Use set_params only when params are fixed or varying!'

    # *** Make versions of show*() that return strings.

    def show_status(self):
        """Print lists of fixed, stepped, and varying params."""
        f = [param.name for param in self.fixed_params]
        s = [param.name for param in self.stepping_params]
        v = [param.name for param in self.varying_params]
        print 'Fixed params:   ', f
        print 'Stepping params:', s
        print 'Varying params: ', v

    def show(self):
        """Print basic parameter info: name, value, doc."""
        # *** Provide explicit support for "derived" params, or
        # expect the user to override this?  Explicit support as
        # AutoNamed descriptors may be desirable to avoid forcing
        # users to implement __init__, though derived params can
        # likely be introduced in on_use rather than __init__.
        # Maybe just try calling a show_extra() method the user can
        # optionally define in their model class.
        for name in self.params:
            print getattr(self, name).show()

    def _score(self, args):
        """
        The quantity *minimized* for inference, as a function of the values
        of any varying parameters, passed in a sequence.  This is a kind of
        "curry" of the score function defined by whatever predictor gets
        mixed in with an inference subclass, intended for use by optimization
        algorithms that require a function with a sequence argument.
        """
        # *** Note this changes *all* varying parameters, even those whose
        # actual value is not changed.  This could cause extra work in the
        # user's on_use method, if it has overhead that is only done when
        # a subset of params are changed.  Perhaps cache param values and
        # check them before setting.  Alternately, require users to do this
        # in on_use if it's an issue.  That's probably safest; it would also
        # handle other situations (e.g., manual alteration of a subset of
        # param values).
        if self.use_unbounded:
            for param, val in zip(self.varying_params, args):
                param._unbounded_set_value(val)
        else:
            for param, val in zip(self.varying_params, args):
                param.set_value(val)
        ## print args, self.objective()
        if self.extremize == minimize:
            return self.score()
        else:
            return -self.score()

    def do_grid(self, vargrid=False, method=None, tol=None, nlog=None):
        """Evaluate the inference on the grid defined by stepped params.

        If there are no varying params, the score function is
        evaluated on the grid and returned in a SimpleValueGrid.
        The vargrid, method and tol arguments are ignored.

        If there are varying params, the score function will be
        extremized w.r.t. them over the grid and this "profile" will
        be returned as a SimpleValueGrid.  If vargrid=True,
        the optimized values of the varying parameters will also be
        collected and returned as a SimpleVectorGrid.  The method and
        tol arguments get passed to the resulting optimize() calls.
        """
        results = ScalarGrid(self.stepping_params)
        vdim = len(self.varying_params)
        #if vargrid and (vdim == 0):
        #    raise InferenceError, 'vargrid argument invalid with no varying params!'
        if self.varying_params and vargrid:
            vresults = VectorGrid(self.stepping_params, self.varying_params)
        opt_val = None
        nstep = 0
        # *** Put try/except around score/optimize so we can return a partial
        # grid if there is a failure after successful calculations.
        while True:
            nstep += 1
            ind = self.step_nums()
            if nlog and nstep%nlog==0:
                print 'Calculating grid step', nstep, '(', ind, ')...'
            # Either optimize or just evaluate.
            if self.varying_params:
                param_vals, extremum, n = self.optimize(method=method, tol=tol)
                results.values[ind] = extremum
                if opt_val is None:
                    opt_val = extremum
                    opt_ind = ind
                else:
                    if self.extremize == maximize:
                        if extremum > opt_val:
                            opt_val, opt_ind = extremum, ind
                    else:
                        if extremum < opt_val:
                            opt_val, opt_ind = extremum, ind
                if vargrid:
                    print vresults.values[ind], param_vals
                    if vdim == 1:
                        vresults.values[ind] = param_vals[0]
                    else:
                        vresults.values[ind] = param_vals
            else:
                results.values[ind] = self.score()
            # Log results if requested.
            if nlog and nstep%nlog==0:
                if self.varying_params:
                    print '  ->', param_vals, extremum, n
                else:
                    print '  ->', results.values[ind]
            try:
                self.next_step()
            except StopStepping:
                break
        if self.varying_params:
            results.opt_val = opt_val
            results.opt_ind = opt_ind
        if vargrid:
            vresults.opt_ind = opt_ind
            return results, vresults
        else:
            return results

    def optimize(self, method=None, tol=None):
        """
        Optimize with respect to varying params.

        The call returns a 3-tuple (params, extremum, extra):
          params = array of best-fit parameter values
          extremum = the optimized statistic (min or max as appropriate)
          extra = any add'l value or values (in a tuple) returned by
                  the specified minimization method (e.g., # iters)
        """
        if method:
            method = method.lower()
        elif self.min_method:
            method = self.min_method.lower()
        else:
            raise InferenceError, 'No minimization method has been specified!'
        if tol is None:
            if self.min_tol is None:
                raise InferenceError, 'No minimization tolerance has been specified!'
        # Most methods require a vector of starting param values and
        # a scale for each; collect this info here (transformed if necessary).
        p, d = [], []
        if self.varying_params == []:
            raise InferenceError, 'No parameters are set to vary for fitting!'
        if self.use_unbounded:
            for param in self.varying_params:
                # Calculate delta for the unbounded version.
                uv = param._unbounded_get_value()
                v = param.get_value()
                try:  # We may need to try a step in both directions.
                    param.set_value(v+param.delta)
                    uv2 = param._unbounded_get_value()
                    d.append(uv2-uv)
                except ParamRangeError:
                    param.set_value(v-param.delta)
                    uv2 = param._unbounded_get_value()
                    d.append(uv-uv2)
                param.set_value(v)
                p.append(uv)
        else:
            for param in self.varying_params:
                p.append(param.get_value())
                d.append(param.delta)
        p = array(p)
        if method == 'powell':
            # Use a diagonal matrix for Powell's start directions.
            dd = identity(len(d)) * array(d, float)
            # For Powell, extra = # of iters
            params, extremum, extra = powell(self._score, p, drxns=dd,
                    ftol=1.e-3, maxit=300)
            # Repeat from the 1st fit, with smaller scales and shuffling
            # the directions to avoid Powell's too-common collapsing to a
            # subspace.
            dd = 0.2 * identity(len(d)) * array(d, float)
            params, extremum, extra = powell(self._score, params[:], drxns=dd,
                    ftol=1.e-3, maxit=200, shuffle=True)
        else:
            raise InferenceError, 'Invalid minimization method specified!'
        if self.extremize == maximize:  # adjust for sign reversal
            extremum = -extremum
        if self.use_unbounded:  # Convert unbounded values back to normal.
            for i, param in enumerate(self.varying_params):
                params[i] = param.get_value()
        return params, extremum, extra


    # *** Potentially confusing that fit returns different types of
    # output for different param status, e.g., just a scalar for
    # all fixed, but [(params), cost] for all varying.

    def fit(self, vargrid=False, method=None, tol=None, nlog=None):
        """
        Calculate a fit or fits of the signal model parameters.

        If no parameters are stepping or varying, just calculate and report
        the value of the score function (fit statistic).

        If a subset of parameters are varying and the rest are fixed, optimize
        with respect to the varying parameters.  Return (params, statistic),
        where params is an array with the optimized parameter values, and
        statistic is the optimized score.

        If there are some stepping parameters and the rest are fixed, evaluate
        the fit on the resulting parameter grid.  Return a SimpleValueGrid
        containing the values of the fit statistic.

        If there are both stepping and varying parameters, optimize
        with respect to the varying parameters at each point on the
        grid (i.e., find the profile fit).  By default, return
        a SimpleValueGrid with the profile fit.  If vargrid is not
        False, return (valuegrid, vargrid), where valuegrid is the
        profile fit SimpleValueGrid, and vargrid is a SimpleVectorGrid
        containing the optimized parameter values at each location
        in the grid.  (Otherwise, vargrid is ignored.)
        """
        if self.varying_params == [] and self.stepping_params == []:
            return self.score()
        elif self.stepping_params == []:
            params, score, iters = self.optimize(method=method, tol=tol)
            self.iters = iters
            return params, score
        else:
            return self.do_grid(vargrid=vargrid, method=method, tol=tol, nlog=nlog)
