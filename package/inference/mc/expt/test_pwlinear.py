from pylab import *
from numpy import *
from scipy.stats import norm
from inference.montecarlo import PWLinear

lo, hi = -10, 10
xvals = linspace(lo, hi, 200)
norm1 = norm(-5., 1.)
norm2 = norm(4., 2.)
yvals = .75*norm1.pdf(xvals) + 1.*norm2.pdf(xvals)

pwl = PWLinear(lo, hi, yvals, cut=0.25)

ion()
plot(xvals, yvals)
plot(xvals, pwl.pdf_vals, 'g:')
plot(pwl.nz_centers, pwl.nz_wts/pwl.nz_widths, 'r.')


samps = pwl.sample(200)
# samps = []
# for i in range(200):
#     samps.append(pwl.sample())
plot(samps, .05*ones_like(samps), 'k.')

xpdf = linspace(-7,7,1000)
plot(xpdf, pwl.pdf(xpdf), 'c-')
show()

# Now test logvals=True.
pwl = PWLinear(lo, hi, log(yvals), logvals=True, cut=0.25)
plot(xvals, yvals)
plot(xvals, pwl.pdf_vals, 'g:')
plot(pwl.nz_centers, pwl.nz_wts/pwl.nz_widths, 'r.')
samps = pwl.sample(200)
plot(samps, .05*ones_like(samps), 'k.')
plot(xpdf, pwl.pdf(xpdf), 'c-')
