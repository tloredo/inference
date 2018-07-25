#!/usr/bin/env python

# To install the package to a directory other than site-packages, use:
#
#     python setup.py install --install-lib <path>


def configuration(parent_package='',top_path=''):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('inference', parent_package, top_path)

    # config.add_subpackage('pie')
    # config.add_subpackage('gauss')
    config.add_subpackage('count')
    # config.add_subpackage('montecarlo')
    # config.add_subpackage('drxn')
    config.add_subpackage('grid')
    # config.add_subpackage('integrate')
    # config.add_subpackage('deriv')
    # config.add_subpackage('signals')
    config.add_subpackage('utils')

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(
        # name = 'inference',  # use if name not given in configuration
        version='0.3',
        maintainer="Tom Loredo",
        maintainer_email="loredo@astro.cornell.edu",
        description="Statistical inference package for Python",
        url="http://www.notyetset.org/",
        license='mixed',
        configuration=configuration)
