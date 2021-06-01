"""
utils:  Utilities for the inference package
"""

from . import numutils
from .numutils import *
from . import ioutils
from .ioutils import *
from . import rng
from .rng import *
from .pl import pl
from .plot_bivar import plot_bivar

# TODO:  Remove dependence on nrmin (in PIE).
from . import nrmin

__version__ = '0.1'

__all__ = ['numutils', 'ioutils', 'rng', 'pl', 'plot_bivar', 'nrmin']
