# This file is executed by __init__.py and ppimport hooks.
"""
Adaptive quadrature algorithms for the Inference package.
"""

from . import adapt, adapt_vec, cub2d

__all__ = ['adapt', 'adapt_vec', 'cub2d']

__doc_title__ = "integrate: Cubature methods for the Inference package."

standalone = True
