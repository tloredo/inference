#!/usr/bin/env python

# To install the package to a directory other than site-packages, use:
#
#     python setup.py install --install-lib <path>


def configuration(parent_package='', top_path=None):
    # we know numpy is a valid import now

    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)
    config.add_subpackage('inference')

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(
        name='inference',
        version='0.3',
        maintainer="Tom Loredo",
        maintainer_email="loredo@astro.cornell.edu",
        description="Statistical inference package for Python",
        url="http://www.notyetset.org/",
        license='mixed',
        packages=['inference'],
        configuration=configuration)
