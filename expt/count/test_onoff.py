from pylab import *
from numpy import *
from inference.count.on_off import OnOff
from inference.grid.hpd import HPD1D
from inference.grid.hpdgrid import HPDGrid1D

def text_phys(x, y, s, fontdict=None, **kwargs):
    """
    A convenience function for rendering text on the current plot axes
    using physical coordinates rather than the default logical (data)
    coordinates.  It left-justifies text horizontally and places the bottom
    at y by default.
    
    This can be accomplished with pylab's :text: function using the
    keyword argument transform=ax.transAxes, where :ax: is the current
    axes object (e.g., from subplot or gca).  text_phys just saves the user
    from having to pass this argument.
    """
    # Add kwargs as required.
    if not kwargs.has_key('horizontalalignment'):
        kwargs['horizontalalignment'] = 'left'
    if not kwargs.has_key('verticalalignment'):
        kwargs['verticalalignment'] = 'bottom'    
    kwargs['transform'] = gca().transAxes
    text(x, y, s, fontdict, **kwargs)

# Parameters for an on-off example.
n_on, T_on = 15, 1.
n_off, T_off = 5, 1.
oo = OnOff(n_on, T_on, n_off, T_off)
s_max = 100.  # upper limit for flat prior on signal

figure(figsize=(12,8))
rc('text', usetex=True)  # Use TeX to render text; needed for some math symbols

#-------------------------------------------------------------------------------
# A "bare" subplot pane for displaying textual info:
subplot(221, frameon=False)
xticks([])
yticks([])
text_phys(.1, .8, r'$N_{\rm on} = %i,   T_{\rm on} = %4.1f$' % (n_on, T_on), fontsize=14)
text_phys(.1, .7, r'$N_{\rm off} = %i,   T_{\rm off} = %4.1f$' % (n_off, T_off), fontsize=14)

# Do a crude model comparison, using a flat prior for the signal strength up
# to s_max.  There isn't much posterior density above s=20, so we just integrate
# up to there.

# Get log (marginal) likelihood values & adjust for the prior factor.
svals = linspace(0, 20, 101)
log_pdf = oo.siglml(svals) - log(s_max)  

# Make an HPD instance.  We'll use it here just for the normalization constant,
# i.e., the marginal likelihood.
hpd1 = HPD1D(log_pdf, 0, 20)
log_marg = hpd1.norm()

# The Bayes factor favoring a nonzero signal is the ratio of the marginal
# likelihood to the likelihood for s=0.
bfac = exp(log_marg - oo.siglml(0.))
text_phys(.1, .5, r'$\hbox{Prior}\; s_{\rm max} = %4.1f \rightarrow B = %4.1f$'\
        % (s_max, bfac), fontsize=14)



#-------------------------------------------------------------------------------
# Plot the ratio of the marginal likelihood to the posterior;
# this checks consistency of their separate algorithms.  The likelihood
# algorithm is approximate by default, though it will end up being exact for
# small numbers of counts; the approximation is useful primarily when the
# on-source counts are large.

subplot(222)
svals = linspace(0, 10, 21)
lml = oo.siglml(svals)  # log marginal likelihood for signal
mp = oo.sigmp(svals)
ratio = exp(lml)/mp
plot(svals, ratio, 'b-')
ylim(0,1.2*mean(ratio))
xlabel('$s$')
ylabel('Ratio')
#text_phys(.5, .85, 'Ratio of likelihood to posterior\n(should be const.)', fontsize=9)
title('Likelihood/posterior ratio')

#-------------------------------------------------------------------------------
# Plot the signal pdf; dashed curve has a longer background integration
# and thus a more precise background estimate.

subplot(223)
fac = 10.  # 2nd case has background integration fac * default case
svals = linspace(0, 20, 101)  # signal levels for plot

# Plot normalized pdf for 1st case.
pdf1 = oo.sigmp(svals)
plot(svals, pdf1, 'b-', label=r'$T_{\rm off}=%4.1f$' % T_off)
# Find 95% HPD region and display it as a transparent fill.
hpd1 = HPD1D(log(pdf1), 0, 20)  # uses already calculated grid; note *log*;
                                # use HPDGrid1D if grid not already available
level, range = hpd1.critlevel(.95)  # returns level & array of boundaries
# For this unimodal case, we expect 1 boundary if the region includes s=0,
# and two boundaries otherwise; for plotting, we need to add the s=0 boundary
# if needed.
if len(range)==1:
    range = (0., range[0])
sfill = linspace(range[0], range[1], 50)
xf, yf = poly_between(sfill, 0., oo.sigmp(sfill))  # so fill extends to y=0
fill(xf, yf, 'b', alpha=0.3)

# Plot the pdf for 2nd case, and fill its 95% HPD region.
oo2 = OnOff(n_on, T_on, int(fac*n_off), fac*T_off)
plot(svals, oo2.sigmp(svals), 'g:', label=r'$T_{\rm off}=%4.1f$' % (fac*T_off,))
# Here we use HPDGrid1D, just to illustrate.
def logpdf2(s):
    return log(oo2.sigmp(s))
hpd2 = HPDGrid1D(logpdf2, 0, 20, 101)
level, range = hpd2.critlevel(.95)
if len(range)==1:
    range = (0., range[0])
sfill = linspace(range[0], range[1], 50)
xf, yf = poly_between(sfill, 0., oo2.sigmp(sfill))
fill(xf, yf, 'g', alpha=0.3)

xlabel('$s$')
ylabel('$p(s\mid D)$')
legend()
#text_phys(.5, .85, 'Marginal pdfs for signal', fontsize=9)
title(r'Marginal pdfs for $s$ \& 95\% HPD regions')

#-------------------------------------------------------------------------------
savefig('test_onoff.png')
show()
