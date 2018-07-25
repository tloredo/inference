from numpy import zeros, empty, linspace, nan
from param import *
from inference.utils import pl, plot_bivar

class ParamValues(object):
    """Container for copies of parameter values, accessible as attributes."""

    # Since a parameter value may be mutable (e.g., an array) and
    # we don't want subsequent changes in the original parameter to affect
    # values in an earlier ParamValues instance, ParamValues keeps
    # simple copies of stored parameter values.  (The copying is
    # implemented by the parameter's copy_value() method.)

    # *** Is there a use case that would require a deep copy?

    # *** Perhaps don't allow altering of attributes.

    def __init__(self):
        self.names = []  # ??? Should this be 'params' for consistency?
        self.docs = {}

    def store(self, name, param):
        if name in self.names:
            raise ParamError, 'Parameter value already stored!'
        self.names.append(name)
        self.docs[name] = param.doc
        setattr(self, name, param.copy_value())

    def __str__(self):
        s = ''
        for name in self.names:
            s += name + " = " + str(getattr(self,name))
            doc = self.docs[name]
            if doc:
                s += ' (' + doc +')'
            s += '\n'
        return s


# Simple grids for storing inference results.
# "Nonsimple" would be "multigrids" that allow one to cover the
# parameter space with patches of various size and resolution,
# like the 1-d self-refining grid we use for frequency spectrum analysis.
# *** Multigrids are not yet implemented.


class ScalarGrid(object):
    """Container for scalar results from computation on stepped param grids.
    E.g., use this for storing values of chi**2, the log likelihood, or the log
    posterior density."""

    def __init__(self, stepping_params):
        self.ndim = len(stepping_params)
        self.axes = []
        self.shape = []
        for param in stepping_params:
            self.axes.append( (param.name, param.slo, param.shi, \
                               param.steps, param.stype) )
            self.shape.append(param.steps)
        self.shape = tuple(self.shape)
        self.values = empty(self.shape, float)
        # Fill values with NaN in case it gets only partly filled due to a
        # failure during a grid calculation.
        self.values[:] = nan

    def __str__(self):
        """Dummy string method."""
        s = str(self.ndim) + '-D ScalarGrid instance:\n'
        for name,lo,hi,steps,stype in self.axes:
            if stype == lin_steps:
                s += '  %s in [%g,%g]; %d linear steps\n' % (name, lo, hi, steps)
            elif stype == log_steps:
                s += '  %s in [%g,%g]; %d log steps\n' % (name, lo, hi, steps)
        return s

    # *** Define __getattr__ to access axes values by name.
    # And/or provide a dictionary of axes and 'name_vals' vectors?

    def plot(self, gscale=False, *args, **kwds):
        """
        Implement default plots for 1-D and 2-D grids of scalar values.
        """
        if self.ndim == 1:
            print 'Plotting 1-D grid values.'
            name, lo, hi, steps, stype = self.axes[0]
            if stype == lin_steps:
                xvals = linspace(lo, hi, steps)
                try:
                    kwds['overlay']  # don't start new figure
                except:
                    pl.figure()  # start a new figure by default
                pl.plot(xvals, self.values, 'b-')
                pl.xlabel(name)
                try:
                    ylab = kwds['ylabel']
                except KeyError:
                    ylab = 'Log Score'
                pl.ylabel(ylab)
            else:
                raise NotImplementedError('Plotting log steps not yet implemented!')
        elif self.ndim == 2:
            print 'Plotting contours of 2-D grid values.'
            x_name, x_lo, x_hi, x_steps, x_type = self.axes[0]
            y_name, y_lo, y_hi, y_steps, y_type = self.axes[1]
            if x_type == lin_steps and y_type == lin_steps:
                xvals = linspace(x_lo, x_hi, x_steps)
                yvals = linspace(y_lo, y_hi, y_steps)
                try:
                    kwds['overlay']  # don't start new figure
                    del kwds['overlay']  # don't pass this on to plot_bivar
                except:
                    pl.figure()  # start a new figure by default
                plot_bivar(xvals, yvals, self.values, xlabel=x_name,
                           ylabel=y_name, **kwds)
            else:
                raise NotImplementedError('Plotting log steps not yet implemented!')
        else:
            raise ValueError('No default plot for grids with dimen > 2!')


class VectorGrid(object):
    """Container for array-valued results from computation on stepped param
    grids, e.g., a set of parameter estimates (of varying_params) throughout the
    grid."""

    # *** Implementation unfinished!

    def __init__(self, stepping_params, stored_params):
        self.ndim = len(stepping_params)
        self.vdim = len(stored_params)
        self.axes = []
        self.shape = []
        for param in stepping_params:
            self.axes.append( (param.name, param.slo, param.shi, \
                               param.steps, param.stype) )
            self.shape.append(param.steps)
        if self.vdim > 1:
            self.shape.append(self.vdim)
        self.shape = tuple(self.shape)
        self.values = zeros(self.shape, float)

    def __str__(self):
        """Dummy string method."""
        s = str(self.ndim) + '-D SimpleVectorGrid instance:\n'
        for name,slo,shi,steps,stype in self.axes:
            if stype == lin_steps:
                s += '  %s in [%g,%g]; %d linear steps\n' % (name, slo, shi, steps)
            elif stype == log_steps:
                s += '  %s in [%g,%g]; %d log steps\n' % (name, slo, shi, steps)
        return s

