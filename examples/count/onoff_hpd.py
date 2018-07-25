from matplotlib.pyplot import *
from numpy import *
from inference.count.onoff import OnOff
from inference.grid.hpd import HPD1D
from inference.grid.hpdgrid import HPDGrid1D


ion()


def text_rel(x, y, s, fontdict=None, **kwargs):
    """
    A convenience function for rendering text on the current plot axes
    using relative (0:1) coordinates rather than the default logical (data)
    coordinates.  It left-justifies text horizontally and places the bottom
    at y by default.

    This can be accomplished with pylab's :text: function using the
    keyword argument transform=ax.transAxes, where :ax: is the current
    axes instance (e.g., from subplot or gca).  text_rel just saves the user
    from having to pass this argument, which can be a nuisance because of
    the need for the axes instance.
    """
    # Add kwargs as required.
    if 'horizontalalignment' not in kwargs:
        kwargs['horizontalalignment'] = 'left'
    if 'verticalalignment' not in kwargs:
        kwargs['verticalalignment'] = 'bottom'
    kwargs['transform'] = gca().transAxes
    text(x, y, s, fontdict, **kwargs)


# Parameters for an on-off example:
n_on, T_on = 15, 1.
n_off, T_off = 5, 1.
oo = OnOff(n_on, T_on, n_off, T_off)

# Upper limit for flat prior on the signal rate:
s_max = 100.

# This is the upper limit for the grid we'll calculate L(s) on; it is
# smaller than the prior limit, s_max, because the likelihood becomes
# negligible for these data well before s_max; no sense wasting effort there.
s_hi = 30.

# We'll compare with a 2nd case that has a background integration interval
# bgfac * default case; we'll give it bgfac*n_off background counts.
bgfac = 10.
oo2 = OnOff(n_on, T_on, int(bgfac*n_off), bgfac*T_off)

# Some overall figure settings:
fig = figure(figsize=(12,5))  # wide figure, for two panels
rc('text', usetex=True)  # use TeX to render text; needed for some math symbols


# A "bare" subplot pane for displaying textual info:
fig.add_axes([0., 0., .35, 1.], frameon=False)  # left, bottom, width, height
xticks([])
yticks([])

# Write the on/off data.
text_rel(.15, .8, r'$N_{\rm on} = %i,   T_{\rm on} = %4.1f$' % (n_on, T_on), fontsize=14)
text_rel(.15, .75, r'$N_{\rm off} = %i,   T_{\rm off} = %4.1f$' % (n_off, T_off), fontsize=14)
text_rel(.15, .65, r'2nd case has $(N_{\rm off},T_{\rm off})\times$%-4.0f' % bgfac, fontsize=14)


# Do a crude model comparison, using a flat prior for the signal strength up
# to s_max.

# First, get log (marginal) likelihood values & adjust for the prior factor.
svals = linspace(0, s_hi, 101)
log_pdf = oo.siglml(svals) - log(s_max)

# Make an HPD instance.  We'll use it here just for the normalization constant,
# i.e., the marginal likelihood.  We could also just use something from
# scipy.integrate.
hpd1 = HPD1D(log_pdf, 0, s_hi)
log_marg = hpd1.norm()

# The Bayes factor favoring a nonzero signal is the ratio of the marginal
# likelihood to the likelihood for s=0.
bfac = exp(log_marg - oo.siglml(0.))

# Do the same for the 2nd case.
log_pdf = oo2.siglml(svals) - log(s_max)
hpd1 = HPD1D(log_pdf, 0, s_hi)
log_marg = hpd1.norm()
bfac2 = exp(log_marg - oo2.siglml(0.))


line = r'$\hbox{Prior}\; s_{\rm max} = %4.1f \rightarrow B = %4.1f, %4.1f$' %\
       (s_max, bfac, bfac2)
text_rel(.15, .55, line, fontsize=14)


# Plot the signal pdf; dashed curve has a longer background integration
# and thus a more precise background estimate.

fig.add_axes([.4, .1, .55, .8])
svals = linspace(0, s_hi, 101)  # signal levels for plot

# Plot normalized pdf for 1st case.
pdf1 = oo.sigmp(svals)
plot(svals, pdf1, 'b-', label=r'$T_{\rm off}=%4.1f$' % T_off)

# Find 95% HPD region and display it as a transparent fill.
hpd1 = HPD1D(log(pdf1), 0, s_hi)  # uses already calculated grid; note *log*;
                                # use HPDGrid1D if grid not already available
level, range = hpd1.critlevel(.95)  # returns level & array of boundaries

# For this unimodal case, we expect 1 boundary if the region includes s=0,
# and two boundaries otherwise; for plotting, we need to add the s=0 boundary
# if needed.
if len(range) == 1:
    range = (0., range[0])
sfill = linspace(range[0], range[1], 50)
fill_between(sfill, 0., oo.sigmp(sfill), color='b', alpha=0.3)

# Plot the pdf for 2nd case, and fill its 95% HPD region.
plot(svals, oo2.sigmp(svals), 'g:', lw=1.5,
     label=r'$T_{\rm off}=%4.1f$' % (bgfac*T_off,))


# Here we use HPDGrid1D, just to illustrate.

def logpdf2(s):
    return log(oo2.sigmp(s))


hpd2 = HPDGrid1D(logpdf2, 0, s_hi, 101)
level, range = hpd2.critlevel(.95)
if len(range) == 1:
    range = (0., range[0])
sfill = linspace(range[0], range[1], 50)
fill_between(sfill, 0., oo2.sigmp(sfill), color='g', alpha=0.3)

xlabel('$s$', fontsize=16)
ylabel('$p(s\mid D)$', fontsize=16)
legend()
title(r'Marginal pdfs for $s$ with 95\% HPD regions')


savefig('onoff_hpd.png')
