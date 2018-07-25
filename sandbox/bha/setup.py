#!/usr/bin/env python

def configuration(parent_package='',top_path=None):
	"""Return the configuration for the "bha" package.  It will
	contain a single module, the C extension built from bhacalc.c."""

	# Create a package configuration instance named "bha".
	from numpy.distutils.misc_util import Configuration
	config = Configuration('bha', parent_package, top_path)

	# Add a C extension to the package, named "bhacalc", and
	# comprised by a single C file, "bhacalc.c".
	config.add_extension('bhacalc', ['bhacalc.c'])

	return config

if __name__ == "__main__":
	from numpy.distutils.core import setup
	setup(**configuration(top_path='').todict())
