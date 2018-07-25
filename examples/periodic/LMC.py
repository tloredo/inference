"""
Analyze arrival time data from the LCM X-ray pulsar, using 
piecewise-constant (PC) light curve models (frequentist & Bayesian)
and the Rayleigh statistic.  A more complete PC analysis of these
data appears in Gregory & Loredo (1996).

The data are X-ray photon arrival times (barycenter-corrected) from
ROSAT observations of PSR 0540-693

BAYESIAN PERIODIC SIGNAL DETECTION: ANALYSIS OF ROSAT OBSERVATIONS OF
PSR 0540-693
Phil Gregory & Tom Loredo
THE ASTROPHYSICAL JOURNAL, 473:1059-1066, 1996 December 20
"""

from pylab import *
from numpy import *
from periodic import *

# Load the LMC pulsar time series from a previously stored pickle.
# This is just an array of arrival times.  See "read_lmc.py".
# lmc = load('lmc.pkl')  # Note we use numpy's load; pylab has one, too.

# Alternately, just load it directly from the data file (slower).
lmc = fromfile('lmc-times.dat',Float,sep=' ')
print 'Loaded', len(lmc), 'event times spanning', lmc[-1]-lmc[0],'s; here we go!'

# Known values of f, fdot --- we'll do a targeted period search (quicker!).
f0 = 19.8528780560
fdot = -1.90010000E-10

# Change to angular frequency.
w0, wdot = 2*pi*f0, 2*pi*fdot

# Make an object supporting piecewise-constant models.
# It supports both conventional chi**2 "epoch folding" and
# the Gregory-Loredo Bayesian model.
pcm = PCModel(lmc, qpts=100)  # Use 100 quadrature points for phase marginalizations

# Define a frequency grid.
nw = 100
dw = linspace(-2.5e-5, 2.5e-5, nw)
wvals = dw+w0

# Do the periodogram calculations.
pcm.set_nbins(2)  # Two-bin lightcurve
lmlvals, chivals = zeros(nw,Float), zeros(nw,Float)
for i, wval in enumerate(wvals):
    if (i+1)%5 == 0: print 'Freq #', i+1, 'of', nw
    lmlvals[i] = pcm.sml_w(wval, wdot) # Scaled marg. likelihood for w
    chivals[i] = pcm.avgchi             # Phase-averaged chi**2

# For plotting on the same axes, scale lml to 2*chi**2.
chimax2 = 2*max(chivals)
p = chimax2 * lmlvals/lmlvals.max()

# Plot both the log posterior and 2*chi**2 (like a log likelihood).
plot(dw, p)
plot(dw,2*chivals)
xlabel(r'$\delta\omega$')
ylabel(r'$\rm{log} p, 2\chi^2$')
show()

# Calculate and plot the Rayleigh statistic on a wider grid.
# Note it (currently) doesn't handle a frequency derivative.
wvals, Rvals = rsgrid(1000, w0-2.5e-4, w0+2.5e-4, lmc)
plot(wvals-w0, Rvals, 'r-')
show()
