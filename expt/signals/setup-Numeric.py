#!/usr/bin/env python
# To build the extension (shared library will be in a "lib..."
# subdirectory of the "build" directory, created if not already present):
#     python setup.py build
# To install the extension in your site-packages directory:
#     python setup.py install
# The latter will build it first if it hasn't been built already.
# You may delete the build directory after installation if it was
# created by this module.

import distutils, os
from distutils.core import setup, Extension

setup(name = "kepler_vec",
    version = "0.1",
    maintainer = "Tom Loredo",
    maintainer_email = "loredo@astro.cornell.edu",
    description = "Kepler radial velocity function",
    url = "http://www.astro.cornell.edu/staff/loredo/python/",
	py_modules = ['kepler'],
    ext_modules = [
    	Extension('fkepler', ['fkeplermodule.c'])
         ]
    )

