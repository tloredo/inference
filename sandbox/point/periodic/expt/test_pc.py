from scipy import *
from numpy import logaddexp
import pylab as pl
from periodic import PCModel
import _pcperiodic as gl

lmc = load('lmc.pkl')
f0 = 19.8528780560  # Best freq for 2 bin model from pcfreqd
fdot = -1.90010000E-10
w0, wdot = 2*pi*f0, 2*pi*fdot

pcm = PCModel(lmc)
pcm.set_nbins(2)
pcm.set_qpts(100)

print pcm.lml_pw(0.,w0,wdot)
print pcm.lml_pw(3.,w0,wdot)
print pcm.lml_pw(0.,w0+1.e-3,wdot)

def pcd(dw, qpts=50):
    """Evaluate sml vs. offset from best freq."""
    w = w0+dw
    pcm.set_qpts(qpts)
    sml = pcm.sml_w(w)
    avgchi = pcm.avgchi
    pcm.set_qpts(0)
    sml2 = pcm.sml_w(w)
    print sml, log(sml) - pcm.offset, avgchi
    print sml2, log(sml2) - pcm.offset, pcm.avgchi

def pcd1(dw, qpts=50):
    w = w0+dw
    pcm.set_qpts(qpts)
    sml = pcm.sml_w(w)
    avgchi = pcm.avgchi
    print sml, log(sml) - pcm.offset, avgchi

def logtrapz(logy, x=None, dx=1.0):
    """
    Integrate y vs x via the trapezoid rule, with log(y) inputs and
    the log of the quadrature as output.
    """
    n_intvls = logy.shape[0]-1
    loghalf = log(.5)
    if x is not None:
        logdel = x[1:] - x[0:-1]
    else:
        logdel = ones(n_intvls)*dx
    logdel = log(logdel)
    lo = logy[0] + loghalf + logdel[0]
    hi = logy[-1] + loghalf + logdel[-1]
    lsum = logaddexp(lo, hi)
    for i in xrange(1,n_intvls):
        lsum = logaddexp(lsum, logy[i] + logdel[i])
    return lsum


# Test logtrapz:
if 0:
    xvals = linspace(1., 3., 100)
    yvals = log(array([x**2 for x in xvals]))
    print exp(logtrapz(yvals, xvals)), 3**3/3. - 1/3.
    import sys
    sys.exit()

pcd(0)

phases = linspace(0,2*pi,51)
chis = []
for phase in phases:
    chis.append(pcm.chi_pw(phase,w0))
chis = array(chis)

wlop, whip = w0-2*pi*3.e-4, w0+2*pi*3.e-4  # prior freq range is 6e-4 Hz
wlo, whi = -2.e-5, 2.e-5
nw = 11
dw = linspace(wlo, whi, nw)
wvals = dw+w0
smlvals, chivals = [], []
for i, wval in enumerate(wvals):
    smlvals.append(pcm.sml_w(wval, wdot))
    chivals.append(pcm.avgchi)
    if (i+1)%5 == 0: print 'Freq # %i of %i...' % (i+1, nw)
smlvals = array(smlvals)
# Estimate marginal likelihood for full freq range.
wvals = w0 + dw
log_marg_like = log(trapz(smlvals/wvals, wvals) - log(whip/wlop)) - pcm.offset
print 'log marginal likelihood:', log_marg_like
chimax = max(chivals)
p = chimax * smlvals/smlvals.max()

pl.plot(dw,p)
pl.plot(dw,chivals)
pl.plot([dw[0],dw[-1]], [chimax,chimax])
pl.xlim(wlo, whi)
pl.show()
