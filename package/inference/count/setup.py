#!/usr/bin/env python

# To install the package to a directory other than site-packages, use:
#
#     python setup.py install --install-platlib <path>

from sys import platform


if platform == "darwin":
    # NumPy's f2py needs these on macOS:
    extra_link_args = ['-undefined', 'dynamic_lookup', '-bundle']
else:
    extra_link_args = None


def configuration(parent_package='',top_path=None):

    # Create a package configuration instance named "example".
    from numpy.distutils.misc_util import Configuration
    config = Configuration('count', parent_package, top_path)

    config.add_extension('_cbmlike', ['_cbmlike.pyf', '_cbmlike.f'],
                         extra_link_args=extra_link_args)

    return config


# The following "boilerplate" actually executes the setup command,
# using as its input a Python dictionary created by the configuration
# object created by the function above.
if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(version='0.3',
          description='count: Tools for inference with counting processes',
          author='Tom Loredo',
          author_email='loredo@astro.cornell.edu',
          license='SciPy License (BSD Style)',
          # py_modules = ['adapt'],
          **configuration(top_path='').todict())
