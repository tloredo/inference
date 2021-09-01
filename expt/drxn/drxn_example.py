"""
Simple examples of Fisher directions --- directions with an
associated axially symmetric error distribution.
"""

from numpy import *
from drxn import *

# These are 27.4 deg apart; expect odds ~ 1 for sig=10 deg.
d1 = FDirection(10.,  0., 90., deg=True)
d2 = FDirection(10.,  0., 62.6, deg=True)

# These are 35.7 deg apart; expect odds ~ 1 for sig=20 deg.
d3 = FDirection(20.,  0., 90., deg=True)
d4 = FDirection(20.,  0., 44.8, deg=True)

print 'Bayes factors for B~1 Gaussian limits, 10deg, 20deg:'
print exp(doublet_lbf(d1,d2)), exp(doublet_lbf(d3,d4))
print

print 'Bayes factors for coincident directions, sig=10 & 20 deg:'
print exp(doublet_lbf(d1,d1)), exp(doublet_lbf(d3,d3))
print

print 'Comparing doublet and n=2 multiplet:'
print doublet_lbf(d1,d2), multiplet_lbf((d1,d2))
print

print 'Comparing triplet and n=3 multiplet:'
print triplet_lbf(d1,d2,d3), multiplet_lbf((d1,d2,d3))

# These are 26 deg apart; p-value of 5%
d5 = FDirection(10.,  0., 90., deg=True)
d6 = FDirection(10.,  0., 64, deg=True)
d7 = FDirection(25.,  0., 90., deg=True)
d8 = FDirection(25.,  0., 64, deg=True)

print 'odds for 26 deg sepn (these have 5% sign. level):'
print 'sig=10:', exp(doublet_lbf(d5,d6))
print 'sig=25:', exp(doublet_lbf(d7,d8))
