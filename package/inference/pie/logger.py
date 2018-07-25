"""
logger:  A stream logger for PIE, logging to stderr.
"""

import logging

__all__ = ['log_debug', 'log_warning']

pielog = logging.getLogger('PIE')
handler = logging.StreamHandler()
# formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
formatter = logging.Formatter('+++ %(levelname)s %(message)s')
handler.setFormatter(formatter)
pielog.addHandler(handler)
pielog.setLevel(logging.WARNING)

def log_debug():
    pielog.setLevel(logging.DEBUG)

def log_warning():
    pielog.setLevel(logging.WARNING)

