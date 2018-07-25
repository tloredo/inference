#!/usr/bin/env python

# To install the package to a directory other than site-packages, use:
#
#     python setup.py install --install-platlib <path>
#
# To install the package in-place (e.g., for debugging):
#
#     python setup.py build
#
# Then copy any built libraries from the build/lib... directory to the working
# directory.

def configuration(parent_package='',top_path=None):

	# Create a package configuration instance named "example".
	from numpy.distutils.misc_util import Configuration
	config = Configuration('periodic', parent_package, top_path)

	config.add_extension('_rayleigh', ['_rayleigh.c'])
	config.add_extension('_pcperiodic', ['_pcperiodic.pyf', '_pcperiodic.f'])

	return config

# The following "boilerplate" actually executes the setup command,
# using as its input a Python dictionary created by the configuration
# object created by the function above.
if __name__ == "__main__":
	from numpy.distutils.core import setup
	setup(version='0.1',
          description='periodic:  Periodic point process analysis tools',
          author='Tom Loredo',
          author_email = 'loredo@astro.cornell.edu',
          license = 'SciPy License (BSD Style)',
          # py_modules = ['periodic'],
          **configuration(top_path='').todict())
