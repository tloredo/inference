#!/usr/bin/env python

def configuration(parent_package='',top_path=None):

	# Create a package configuration instance named "fkepler".
	from numpy.distutils.misc_util import Configuration
	config = Configuration('signals', parent_package, top_path)

	# Add a C extension to the package, named "fkepler", and
	# comprised of a single C file, "fkeplermodule.c".
	config.add_extension('fkepler', ['fkeplermodule.c'])

	return config

if __name__ == "__main__":
	from numpy.distutils.core import setup
	setup(version='0.1',
          description='signals:  Modules implementing various types of parameterized signal models.',
          author='Tom Loredo',
          author_email = 'loredo@astro.cornell.edu',
          #maintainer_email = 'scipy-dev@scipy.org',
          license = 'SciPy License (BSD Style)',
          **configuration(top_path='').todict())
