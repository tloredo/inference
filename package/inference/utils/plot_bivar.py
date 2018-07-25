"""
Objects supporting plotting of bivariate functions arising in statistics,
e.g., pdfs, log-likelihoods, and chi**2 functions.

Much of the code is adapted from Matplotlib's contour_demo.py.

Adapted from prior code 10 July 2008 by Tom Loredo
"""

from numpy import *
import scipy.stats as stats
from inference.utils import pl

__all__ = ['plot_bivar']


labels = { 'fontsize' : 14 }

crosshair_p = { 'color' : '0.5', 'linestyle' : ':' , 'linewidth' : '1.5'}

def crosshair(ax, x, y):
    ax.axvline(x, **crosshair_p)
    ax.axhline(y, **crosshair_p)

def normal_plevels(n):
    """
    Return an array of values of the probability within a +- k*sigma region
    centered on the mean of a normal distribution, for k=1 to n.
    """
    c1cdf = stats.chi2(1).cdf
    levels = []
    for i in range(1,n+1):
        levels.append(c1cdf(i**2))
    return array(levels)

normal_plevels5 = normal_plevels(5)
normal_plevels6 = normal_plevels(6)
#array([c1cdf(1), c1cdf(4), c1cdf(9), c1cdf(16), c1cdf(25)])

# Ten colors to cycle contours through.
colors=('blue', 'darkgreen', 'firebrick', 'saddlebrown', 'dimgray',
        'k', 'k', 'k', 'k', 'k')

# chi-squared dist'n with 2 dof, for finding default 2-D contour levels.
chisqr_2 = stats.chi2(2)

def plot_bivar(x, y, logpdf, iorder='xy', levels=None, plevels=None, cross=None,
               xlabel='X', ylabel='Y', lw=1, colors=colors, gscale=False,
               gscale_bar=False, clabel=False):
    """
    Plot contours (and optionally a grayscale image and crosshairs) using a
    2-D grid of values of a function evaluated on linearly spaced grids in
    x and y.  Default contour levels are suitable for functions interpretable
    as a bivariate log-pdf with approximately normal contour levels; the hpd
    and hpdgrid modules provide tools for calculating levels accurately
    normalized over the grid.
    
    `iorder` specifies the index order for the `logpdf` array.  If it is 'xy',
    the `logpdf` array is interpreted so that logpdf[i,j] is the value for
    (x[i], y[j]).  Note that this means that x increases along *columns*
    (different choices of row, i).  This is the reverse interpretation
    of Matplotlib's contour and imshow; we use a transposed view of logpdf
    for the plotting.  If `iorder` is 'yx', logpdf[i,j] is the value for
    (y[i], x[j]).

    If `gscale` is not False, it should be a float specifying a factor used
    to set the lower limit (black level) for a grayscale image.  The factor
    multiplies the lowest value of `levels`.  E.g., gscale=10 will produce
    an image with the noticeable gray values highlighting the "tails" of
    `logpdf` (values below the lowest contour).  For backward compatibility,
    using gscale=True corresponds to gscale=10. 
    """
    # *** Add args to handle chi**2 input; internally convert to log-like.
    # This is based on Matplotlib's contour_demo.py
    xg, yg = meshgrid(x, y)
    xlo, xhi = xg[0,0], xg[0,-1]
    ylo, yhi = yg[0,0], yg[-1,0]
    if iorder=='xy':
        logpdf = logpdf.transpose()
    elif iorder != 'yx':
        raise ValueError('Illegal iorder!')
    logmax = logpdf.max()
    dlogpdf = logpdf - logmax
    if levels is None:
        if plevels is None:
            # Base levels on bivariate normal.
            levels = -0.5 * chisqr_2.ppf(normal_plevels6)
        else:
            levels = -0.5 * chisqr_2.ppf(plevels)
    else:
        # Treat levels as *absolute* levels (not relative to max).
        levels = array(levels) - logmax

    # Plot a grayscale image displaying tails beyond the lowest level.
    if gscale:
        if gscale is True:
            gscale = 10.
        lo, hi = gscale*levels.min(), 0.
        im = pl.imshow(dlogpdf, interpolation='bilinear', origin='lower', aspect='auto',
                    vmin=lo, vmax=hi, cmap=pl.cm.gray, extent=(xlo, xhi, ylo, yhi))

    # Plot contours on top.
    cset = pl.contour(xg, yg, dlogpdf, levels, colors=colors, origin='lower',
        linewidths=lw)

    if clabel:
        pl.clabel(cset, fmt='%1.2f', inline=1, fontsize=10)

    # Make a colorbar for the contour lines.
    #CB = colorbar(cset, shrink=0.6, extend='both')

    #title('Doublet Periodogram')
    pl.xlabel(xlabel)
    pl.ylabel(ylabel)
    #hot()  # Now change the colormap for the contour lines and colorbar
    #flag()

    # We can still add a colorbar for the image, too.
    if gscale and gscale_bar:
        CBI = pl.colorbar(im, orientation='horizontal', fraction=.1, shrink=0.8)

    if cross:
        xc, yc = cross
        xlo, xhi = xlim()
        ylo, yhi = ylim()
        pl.hlines(yc, xlo, xhi, 'teal')
        pl.vlines(xc, ylo, yhi, 'teal')
    # This makes the original colorbar look a bit out of place,
    # so let's improve its position.
    #l,b,w,h = gca().get_position().bounds
    #ll,bb,ww,hh = CB.ax.get_position().bounds
    #CB.ax.set_position([ll, b+0.1*h, ww, h*0.8])

    # Use a coordinate formatter that shows the z value in interactive plots.
    # *** This assumes a linearly spaced grid.
    def format_coord(x, y, xlo=xlo, xhi=xhi, ylo=ylo, yhi=yhi, data=dlogpdf):
        ny, nx = data.shape
        dx = (xhi-xlo)/(nx-1.)
        dy = (yhi-ylo)/(ny-1.)
        col = int((x-xlo)/dx+0.5)
        row = int((y-ylo)/dy+0.5)
        if col>=0 and col<nx and row>=0 and row<ny:
            z = data[row,col]
            return 'x=%1.4f, y=%1.4f, z=%1.4f'%(x, y, z)
        else:
            return 'x=%1.4f, y=%1.4f'%(x, y)

    pl.gca().format_coord = format_coord
