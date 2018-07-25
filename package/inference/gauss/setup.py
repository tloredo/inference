#!/usr/bin/env python

def configuration(parent_package='',top_path=None):
	"""Return the configuration for the "vba" package.  It will
	contain a single module, the C extension built from vbamodule.c."""

	# Create a package configuration instance named "vba".
	from numpy.distutils.misc_util import Configuration
	config = Configuration('gauss', parent_package, top_path)

	# Add a C extension to the package, named "vba", and
	# comprised by a single C file, "vbamodule.c".
	config.add_extension('vba', ['vbamodule.c'])

	return config

if __name__ == "__main__":
	from numpy.distutils.core import setup
	setup(version='0.1',
          description='gauss:  Modules for inference tasks using the Gaussian (normal) distribution',
          author='Tom Loredo',
          author_email = 'loredo@astro.cornell.edu',
          license = 'SciPy License (BSD Style)',
          **configuration(top_path='').todict())
