"""
pl:  Pylab plotting module interface module for the inference package.

If matplotlib is available, its 'pylab' module is returned as "pl".
Otherwise, "pl" points to an object that raises an AttributeError whenever
attempts are made to access an attribute or method, with a warning about
missing pylab.

This allows inference modules to 'import pl' without raising import errors
if the user happens to not have matplotlib.
"""


__all__ = ['pl']


class NoPyplot(object):

    def __getattribute__(self, name):
        raise AttributeError("No plot capability --- Matplotlib's pyplot unavailable!")


try:
    import matplotlib.pyplot as pl
except ImportError:
    pl = NoPyplot()

