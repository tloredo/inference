from numpy import *
from _asurvkm import kmestm, kmdif, quart

# Censoring index (-1 = censored):
gal1_ind = array([0, 0, -1, -1, 0, -1])
# IR abs. magnitudes for galaxies:
gal1_mag = array([28.5, 26.9, 29.7, 28.1, 30.1, 27.6])

# Diff'l KM estimate bin params:
nb, w, m_l = 5, 2., 25.

ierr,sx,vx,smean,err,nu,su,nc,sc =  kmestm(gal1_ind,gal1_mag)
sx = sx[:nu]
su = su[:nu]
sc = sc[:nc]
ntot = gal1_mag.shape[0]
print 'KM CDF:'
print su
print sc
print

# Diff'l:
bin_l,bin_u,diff = kmdif(sx, su, ntot, m_l, w, nb) 
print 'Diff\'l:'
print bin_l
print bin_u
print diff
print

# Quartiles:
q = quart(su, sc)
print 'Quartiles:', q