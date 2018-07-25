#!/usr/bin/env python

def configuration(parent_package='',top_path=None):

	# Create a package configuration instance.
	from numpy.distutils.misc_util import Configuration
	config = Configuration('montecarlo', parent_package, top_path)

	# Add C/Fortran extensions to the package.
	config.add_extension('_ppssampler', ['_ppssampler.c', 
        'sampford.c', 'compact5table.c','randomkit.c'])
	config.add_extension('_mvnt', ['_mvnt.pyf', '_mvnt.f'])

	return config

if __name__ == "__main__":
	from numpy.distutils.core import setup
	setup(version='0.1',
          description='montecarlo:  Modules implementing Monte Carlo sampling methods',
          author='Tom Loredo',
          author_email = 'loredo@astro.cornell.edu',
          license = 'SciPy License (BSD Style)',
          **configuration(top_path='').todict())
