from numpy import *
from KM import KM

# ASURV's gal1 example:

# Censoring index (-1 = censored):
gal1_ind = array([0, 0, -1, -1, 0, -1])
# IR abs. magnitudes for galaxies:
gal1_mag = array([28.5, 26.9, 29.7, 28.1, 30.1, 27.6])

km = KM(all=gal1_mag, ind=gal1_ind)
a, o, ea, eo, ee = km.cdf_pts(25,31)
km.diff(5, 25, 2)
