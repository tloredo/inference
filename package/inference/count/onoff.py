"""
Bayesian inference of the signal rate of a Poisson counting process, where
the total rate is signal plus background.

This is based on analytic marginalization of the background rate, as described
by Loredo (1989, 1992).

By Tom Loredo
"""


from numpy import zeros_like
from ._cbmlike import sigctp, psigc, slml_offset, vslml_offset, slmlike, \
    vslmlike, slmliked

__all__ = ['OnOff',
           'slmlike', 'vslmlike', 'slmliked', 'slml_offset', 'vslml_offset',
           'sigctp', 'psigc']


class OnOff(object):
    """
    Model on-source (signal + background) and off-source (background only)
    event count data with (constant) Poisson counting processes, providing
    various Bayesian inferences (method names in parens):

      * The log marginal likelihood for the signal rate (siglml);
      * The marginal posterior density for the signal rate (sigmp);
      * HPD credible regions for the signal (sighpd);
      * The posterior distribution for the number of events attributable to
        the signal (pnsig).

    These inferences all use flat (improper) priors for the signal and 
    background rates, as needed (the marginal likelihood for the signal does
    not include a prior factor).
    """

    def __init__(self, on_cts, on_intvl, off_cts, off_intvl, cutoff=1.e-5):
        """
        Initialize an OnOff object.

        :Parameters:

          on_cts : int
            Counts observed on source (source + background)
          on_intvl: float
            Interval spanned on source
          off_cts: int
            Counts observed off source (background only)
          off_intvl : float
            Interval spanned off source

        :Keywords:

          cutoff : float
            Cutoff for truncation of the sum of counts assigned to the signal
            in the log marginal likelihood calculation.

            This should be in [0., 1.]; terms in the sum less than cutoff * the
            maximum value will be omitted.  Cutoff=0 produces exact results
            (to machine precision).

            cutoff does not affect log marginal posterior calculations, which
            currently always use the full sum.

        :Notes:

            The __init__ arguments are stored and may be accessed as data 
            members, but they should not be modified after object construction.
        """
        self.on_cts = on_cts
        self.on_intvl = on_intvl
        self.off_cts = off_cts
        self.off_intvl = off_intvl
        self.cutoff = cutoff
        self.offset = slml_offset(on_cts, on_intvl, off_cts, off_intvl)
        self.p_signal = None  # only calculate & cache if needed

    def siglml(self, s):
        """
        Log of the (marginal) likelihood for scalar or vector signal rate(s).

        :Parameters:

          s : float OR float array
            Signal rate(s) (counts per unit interval)

        :Returns:

          siglml : float OR float array
            Natural logarithm of the marginal likelihood for the signal rate(s)

        :Exceptions:

          ValueError
            Raised for underflow/overflow in the series; try adjusting
            cutoff if this is an issue
        """
        try:
            s = float(s)  # raises TypeError for arrays of length != 1
            llike, nt, err = slmlike(s, self.on_cts, self.on_intvl,
                                     self.off_cts, self.off_intvl, self.offset, self.cutoff)
            if err != 0:
                raise ValueError('Underflow/overflow in likelihood calculation!')
            return llike
        except TypeError:
            if len(s.shape) != 1:
                raise ValueError('sigll handles only 1-D arrays!')
            llvals = zeros_like(s)
            for i, sval in enumerate(s):
                llvals[i], nt, err = slmlike(sval, self.on_cts, self.on_intvl,
                                             self.off_cts, self.off_intvl, self.offset, self.cutoff)
                if err != 0:
                    raise ValueError('Underflow/overflow in likelihood calculation!')
            return llvals

#         try:
#             # This will raise AttributeError if s is not a Numpy array.
#             if len(s.shape) != 1:
#                 raise ValueError('sigll handles only 1-D arrays!')
#             llvals = zeros_like(s)
#             for i, sval in enumerate(s):
#                 llvals[i], nt, err = slmlike(sval, self.on_cts, self.on_intvl,
#                     self.off_cts, self.off_intvl, self.offset, self.cutoff)
#                 if err != 0:
#                     raise ValueError('Underflow/overflow in likelihood calculation!')
#             return llvals
#         except AttributeError:
#             llike, nt, err = slmlike(s, self.on_cts, self.on_intvl,
#                        self.off_cts, self.off_intvl, self.offset, self.cutoff)
#             if err != 0:
#                 raise ValueError('Underflow/overflow in likelihood calculation!')
#             return llike

    def sigmp(self, s):
        """
        Marginal posterior density for scalar or vector signal rate(s).

        The calculation is exact (to machine precision), and requires a sum
        with a number of terms equal to the number of on-source counts.  In
        many cases with a large number of counts, many of the terms will be
        negligible.  The marginal likelihood method, sigml, can neglect small
        terms and may be preferable to sigmp in these cases; its results
        will have to be explicitly normalized to produce a posterior pdf.

        :Parameters:

          s : float OR float array
            Signal rate(s) (counts per unit interval)

        :Returns:

          siglmp : float OR float array
            Natural logarithm of the marginal posterior for the signal rate(s)
        """
        # Make sure we've calculate an array of signal count probabilities.
        if self.p_signal is None:
            self.p_signal = sigctp(self.on_cts, self.on_intvl,
                                   self.off_cts, self.off_intvl)
        # Treat array or vector argument cases separately.
        try:
            s = float(s)  # raises TypeError for arrays of length != 1
            return psigc(s, self.on_intvl, self.p_signal)
        except TypeError:
            if len(s.shape) != 1:
                raise ValueError('sigll handles only 1-D arrays!')
            llvals = zeros_like(s)
            for i, sval in enumerate(s):
                llvals[i] = psigc(sval, self.on_intvl, self.p_signal)
            return llvals

#         try:
#             # This will raise AttributeError if s is not a Numpy array.
#             print s, s.shape, type(s)
#             if len(s.shape) != 1:
#                 raise ValueError('sigll handles only 1-D arrays!')
#             llvals = zeros_like(s)
#             for i, sval in enumerate(s):
#                 llvals[i] = psigc(sval, self.on_intvl, self.p_signal)
#             return llvals
#         except AttributeError:
#             return psigc(s, self.on_intvl, self.p_signal)

    def pnsig(self, n=None):
        """
        The probability that n events are from the signal source.

        :Parameters:

          n : int
            A number of events <= on_cts

        :Returns:

          p : float OR float array
            If n is provided, the probability that n of the on_cts are from
            the signal source is returned; with no arguments, an array of
            length (on_cts+1) is returned so that p[n] is the probability
            that n counts are from the signal source.
        """
        # Make sure we've calculate an array of signal count probabilities.
        if self.p_signal is None:
            self.p_signal = sigctp(self.on_cts, self.on_intvl,
                                   self.off_cts, self.off_intvl)
        if n:
            if n < 0 or n > self.on_cts:
                raise ValueError('Number of source counts must be <= on_cts!')
            return self.p_signal[n]
        else:
            return self.p_signal[:]  # return copy so cache can't be altered


# Interim documentation for the underlying Fortran extension module.
_cbmlike_doc = \
    """
This module docstring is copied from the underlying Fortran-77
implementation.  It will eventually be revised to better reflect 
the Python interface.  Note that each function in this module
has a docstring detailing its Python interface; access it in
the Python interpreter via "print <function>.__doc__", or more
simply in IPython via "?<function>".

CBMLike:  Constant Background Marginal Likelihood

Subprograms for calculating the Poisson likelihood for
a signal rate when data are collected with a constant
unknown background rate, estimated by observations
in "off-source" intervals.

NOTE:  These docs (at the start of this file) are not
completely updated wrt the following code; trust each
subprogram's docs over this header!

The CBM likelihood is proportional to equation (12.18) of
"The Promise of Bayesian Inference for Astrophysics"
(in *Statistical Challenges in Modern Astronomy*, ed.
E. Feigelson and G. Babu, Springer-Verlag, 1992,
pages 275-305).  That equation gives the normalized marginal
posterior density for the signal.  These subroutines calculate 
the marginal likelihood for the signal.  This is not
normalized over signal values.  It is also not
normalized over data (on/off count) values; it is
missing data-dependent factors that would make it a
rigorously normalized sampling distribution.

The data that the user must provide are:
    n_on        Integer # of counts seen on-source
    t_on        On-source interval size (eg time for
                rate measurements, time*area for flux
                measurements)
    n_off        Integer # of counts seen off-source
    t_off        Off-source interval size

The user must also set two parameters:
    s                The signal strength for which the likelihood
                is to be calculated
    cut                A parameter specifying the accuracy of
                the calculation, in [0., 1.). 

The calculation is done by summing terms in the sum
given in equation (12.18) mentioned above.  If CUT = 0.,
the sum is done with all (n_on + 1) terms, giving an
exact result (to machine precision).  Otherwise, the
sum starts with the largest term (call it LTERM), and
adds decreasing terms until the new terms fall below
CUT*LTERM.  Terms in the sum are calculated recursively.

Note that the C_i factor in the sum is the probability
that i of the n_on counts are actually from the signal,
not the background.  For large values of n_on or n_off,
many C_i coefficients will be negligibly small, possibly
causing an underflow.  There is no use keeping such terms,
so I strongly recommend using CUT > 0.  This obviously
also saves computational time, though in none of my own
applications has this likelihood calculation demanded
excessive time (with, say, CUT = 1.e-5).

When the number of counts on- or off-source is large,
the value of the likelihood can be very small, possibly
smaller than allowed by Fortran77 double precision.
Thus these subroutines calculate a LOGARITHM of the 
SCALED likelihood.  The user must first calculate a scale
factor using SL_SCL, and then calculate the scaled log
likelihood using SLLIKE.  The scale factor depends only
on the data, not on the values of S or CUT, hence the
separate functions.  

Thus, once the data are known, the user should execute
the following:

    double precision sl_scl
    double precision lscl, t_on, t_off
    integer n_on, n_off

    lscl = sl_scl(n_on, t_on, n_off, t_off)

Then for each value of S for which the likelihood is needed,
execute (for example):

    double precision sllike
    double precision s, cut, sll
    integer nt

    cut = 1.e-5
    sll = sllike(s, n_on, ta_on, n_off, ta_off, lscl, cut, nt)

SLL will be the LOG of the SCALED likelihood.  If you need
the actual value of the log likelihood for some reason,
just add LSCL to SLL.  But since the value of LSCL depends
only on the data and not on the value of S, no common
calculations (parameter estimates, credible regions, Bayes
factors, odds ratios) are afffected by using only the
scaled value.

NT is set by the function and indicates the number of terms
used in the calculation.  NT <= N_ON + 1.  It may be of
interest for monitoring how many terms are necessary at
a particular level of accuracy (value of CUT).

If the user knows beforehand that scaling is not necessary
(e.g., the n_on and n_off values are relatively small),
the small amount of overhead associated with sl_scl can be
eliminated by eliminating the call to sl_scl and simply
setting LSCL=0.

If the derivative of the scaled log likelihood wrt S is also
needed, use the following in place of the last set of
Fortran77 lines:

    double precision s, cut, sll, dsll
    integer nt

    cut = 1.e-5
    call cbmlike (s, n_on, ta_on, n_off, ta_off, lscl, cut, 
                  sll, dsll, nt)

SLL is as above; DSLL is its derivative wrt S.
"""

# Redefine docstrings:

slmlike.__doc__ += """
"Signal Marginal Log LIKElihood"

Return the log likelihood for the signal, calculating with all
terms up to a factor of CUT from the leading term.

offset is the log of a scale factor (see slml_offset, above),
used to prevent underflow/overflow.

On return, nt is the number of terms used; nt <= n_on+1.

exp(slmlike) is NOT normalized w.r.t. s!  It's a log LIKELIHOOD, 
not a posterior density.

It is also NOT normalized w.r.t. n_on, n_off.  It is a
likelihood function, missing several data-dependent
factors that would make it a normalized sampling dist'n.
    
See the cbmlike module docstring for more information on useage.
"""

vslmlike.__doc__ += """
"Vector Signal Marginal Log LIKElihood"

Return the log likelihood for the signal values in nc channels
of on/off data, calculating with all
terms up to a factor of CUT from the leading term.

off is the log of a scale factor (see slml_off, above),
used to prevent underflow/overflow.

exp(slml) is NOT normalized w.r.t. s!  It's a log LIKELIHOOD, 
not a posterior density.

It is also NOT normalized w.r.t. n_on, n_off.  It is a
likelihood function, missing several data-dependent
factors that would make it a normalized sampling dist'n.

See the cbmlike module docstring for more information on useage.
"""

slmliked.__doc__ += """
"Signal Marginal Log LIKElihood and Derivative"

Return the scaled log likelihood for the signal (SLL), 
calculating with all terms up to a factor of CUT from the 
leading term.  Also calculate the derivative of SLL
with respect to s (returned as DSLL).

offset is the log of a scale factor (see slml_off, above),
used to prevent underflow/overflow.

On return, nt is the number of terms used; nt <= n_on+1.

SLL is NOT normalized w.r.t. s!  It's a log LIKELIHOOD, not
a posterior density.

noff=0 gives nan

See the cbmlike module docstring for more information on useage.
"""

slml_offset.__doc__ += """
"Constant Background Marginal Log Likelihood OFFset"

Return the log of a scale factor for use in SLIKE, to
prevent overflow/underflow.

See the cbmlike module docstring for more information on useage.
"""

vslml_offset.__doc__ += """
"Vector Constant Background Marginal Log Likelihood OFFSET"

Return the log of a scale factor for use in SLIKE, to
prevent overflow/underflow, for each channel in a
vector of on/off measurements

See the cbmlike module docstring for more information on useage.
"""

sigctp.__doc__ += """
SIGnal CounT Probabilities

Calculate the c coefficients in the series expression for the
for the marginal likelihood of the signal.  The dimension of c
is n_on + 1, and the first element is c(0).

c(i) is the probability that i of the n_on counts are from
the signal.

The coefficients are calculated recursively, going up and
down from the largest one so roundoff error accumulates only
on the negligible terms.

This assumes a flat prior for the background rate.

See the cbmlike module docstring for more information on useage.
"""

psigc.__doc__ += """
Probability for SIGnal from Coefficients

Calculate the probability density for the signal, s, using the
c coefficients from coeff, above.

See the cbmlike module docstring for more information on useage.
"""
