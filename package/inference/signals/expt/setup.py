#!/usr/bin/env python

def configuration(parent_package='',top_path=None):
	"""Return the configuration for the "fkepler" package.  It will
	contain a single module, the C extension built from fkeplermodule.c."""

	# Create a package configuration instance named "fkepler".
	from numpy.distutils.misc_util import Configuration
	config = Configuration('fkepler', parent_package, top_path)

	# Add a C extension to the package, named "fkepler", and
	# comprised of a single C file, "fkeplermodule.c".
	config.add_extension('fkepler', ['fkeplermodule.c'])

	return config

if __name__ == "__main__":
	from numpy.distutils.core import setup
	setup(version='0.1',
          description='kepler - functions facilitation Keplerian orbit calculations',
          author='Tom Loredo',
          author_email = 'loredo@astro.cornell.edu',
          #maintainer_email = 'scipy-dev@scipy.org',
          license = 'SciPy License (BSD Style)',
          **configuration(top_path='').todict())
