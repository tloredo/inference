from numpy import pi, zeros, log, linspace
from scipy.special import gammaln
import _pcperiodic as gl  # aka "Gregory-Loredo algorithm"
from _rayleigh import rayleigh, rsgrid, lml_wk

__all__ = ['PCModel', 'rayleigh', 'rsgrid', 'lml_wk']

twopi = 2*pi

class PCModel(object):
    """
    Piecewise-constant model for a periodic Poisson point process.
    "Periodic" is understood somewhat loosely to include processes with
    possibly nonzero first and second frequency derivatives.

    In the method names for this class, "w" stands collectively for
    angular frequency (w) and its derivatives (wdot, wddot), and "p" for
    phase (over [0,2*pi]).
    """

    def __init__(self, data, nbins=None, qpts=None):
        """
        Initialize a piecewise-constant model, copying the data, and optionally
        setting the number of bins for the PC model, and
        the number of points for a quadrature marginalizing over phase
        (these can be altered after initialization).  qpts=0 results in an
        exact calculation.
        """
        if len(data.shape) != 1:
            raise ValueError('Input must be 1-d float array!')
        self.data = data.copy()
        self.ndata = len(data)
        self.lndfac = gammaln(self.ndata+1)
        if qpts is None:
            self.qpts = None
            self.absc = None
        elif qpts >= 0:
            self.qpts = qpts
            if nbins and qpts > 0:
                # This depends on nbins b/c we integ. over one bin width.
                self.absc = gl.trap_absc(qpts, nbins)
            else:
                self.absc = None
        else:
            raise ValueError('qpts must be a positive integer!')
        if nbins is not None:
            self.set_nbins(nbins)
        else:
            self.bins = None
        self.offset = None # Offset for logs to avoid over/underflow
        self.w = None
        self.wdot = 0.
        self.wddot = 0.
        self.phases = zeros(len(data))

    def set_nbins(self, nbins):
        """Set the number of bins for the piecewise-constant rate model."""
        self.nbins = nbins
        self.bins = zeros(nbins, int)  # For storing event counts
        # *** Rewrite .f so these two could be replaced by one call.
        self.phibins = zeros(nbins)
        gl.binbounds(self.phibins)  # Phase boundaries for bins
        self.offset = None
        # nbins-dependent part of the Ockham factor:
        self.logbfac = gammaln(nbins) + self.ndata*log(nbins) + \
                       gammaln(self.ndata+1) - gammaln(self.ndata+nbins)
        if self.qpts > 0:
            self.absc = gl.trap_absc(self.qpts, nbins)
        else:
            self.absc = None

    def set_qpts(self, n):
        if n >= 0:
            self.qpts = n
            if self.nbins and n > 0:
                self.absc = gl.trap_absc(n, self.nbins)
            else:
                self.absc = None
        else:
            raise ValueError('qpts must be a positive integer!')

    def _set_offset(self):
        # Bin data with current w params, and find the
        # multiplicity for 0 phase.
        gl.fbindata(self.phases, 0., self.phibins, self.bins)
        self.offset = - gl.limult(self.ndata, self.bins)

    def _update_params(self, w, wdot, wddot):
        """
        Update frequency parameters, and intermediate
        results that are stored for future use, as necessary.
        Unspecified parameters (val=None) retain any previous settings.
        """
        if self.nbins is None:
            raise ValueError('Must specify number of bins for the model!')
        #print 'new:', w, wdot, wddot
        #print 'old:', self.w, self.wdot, self.wddot
        # If params changed, calculate event phases.
        if (w != self.w) or (wdot != self.wdot) or (wddot != self.wddot):
            # Require w be specified, in this call or a previous one.
            if w is None:
                if self.w is None: raise ValueError('Need a w parameter!')
                w = self.w
            else:
                self.w = w
            if wdot is not None:
                self.wdot = wdot
            else:
                wdot = self.wdot
            if wddot is not None:
                self.wddot = wddot
            else:
                wddot = self.wddot
            if wddot is None:
                if wdot is None:
                    gl.phases0(w, self.data, self.phases)
                else:
                    gl.phasesd(w, wdot, self.data, self.phases)
            else:
                gl.phasesdd(w, wdot, wddot, self.data, self.phases)

    def lml_pw(self, phi, w=None, wdot=None, wddot=None):
        """
        Log marginal likelihood for phase and angular frequency parameters,
         marginalized over signal shape.  This is just the log of the
        inverse multiplicity of the counts in bins.

        Add self.logbfac to make it additionally a likelihood over nbins.

        Since this is fast (and seldom needed), we don't do any bookkeeping
        of previous binning that may have used the same parameters.
        """
        self._update_params(w, wdot, wddot)
        # This uses phi and phibins to bin the events.
        gl.fbindata(self.phases, phi, self.phibins, self.bins)
        lml = gl.limult(self.ndata, self.bins)
        if self.offset is None: self.offset = - lml
        return lml

    def chi_pw(self, phi, w=None, wdot=None, wddot=None):
        """
        Chi**2 vs. phase and angular frequency parameters.
        """
        self._update_params(w, wdot, wddot)
        gl.fbindata(self.phases, phi, self.phibins, self.bins)
        return gl.chisqr(self.ndata, self.bins)

    def sml_w(self, w=None, wdot=None, wddot=None):
        """
        Scaled marginal likelihood for angular frequency parameters & bin
        number, marginalized over phase and signal shape.

        The actual marginal likelihood is sml_w*exp(-self.offset).  Since one
        might integrate this over millions of w values, the result is left
        scaled so the user can sum the returned sml_w values.  The log of
        the integral minus self.offset is the log of the full marginal
        likelihood.
        """
        self._update_params(w, wdot, wddot)
        if self.offset is None: self._set_offset()
        if self.qpts==0:
            self.lml, self.avgchi = gl.qlw_exact(self.lndfac,
                self.phases, self.phibins, self.bins, self.offset)
        else:
            self.lml, self.avgchi = gl.qlw_trap(self.lndfac,
                self.phases, self.phibins, self.bins, self.absc, self.offset)
        return self.lml

    def avgchi_w(self, w=None, wdot=None, wddot=None):
        """
        Phase-averaged chi**2 vs. frequency parameters.
        """
        self._update_params(w, wdot, wddot)
        if self.offset is None: self._set_offset()
        if self.qpts==0:
            self.lml, self.avgchi = gl.qlw_exact(self.lndfac,
                self.phases, self.phibins, self.bins, self.offset)
        else:
            self.lml, self.avgchi = gl.qlw_trap(self.lndfac,
                self.phases, self.phibins, self.bins, self.absc, self.offset)
        return self.avgchi

    def lml_const(self):
        """
        Return the log marginal likelihood for the constant model.
        """
        return gl.lodds()

    def rate_est(self, phases, w=None, wdot=None, wddot=None):
        """
        Estimate the event rate as a function of phase through one period,
        for a model with the given frequency parameters and bin number,
        marginalizing over phase.  The phase is in [0,2*pi].

        If phases is an integer, the rate is calculated at that number of
        phases over a period.  A tuple of (phases, rates, sigmas) is
        returned, each member being a 1-d array of values.

        If phases is an array or sequence of numbers, the rate is calculated
        at those phases, and only (rates, sigmas) is returned.
        """
        self._update_params(w, wdot, wddot)
        if self.offset is None: self._set_offset()
        try:
            # phases is an array or sequence
            ret_phases = False
            len(phases)
            phases = array(phases)
        except AttributeError:  # Raised if len(phases) fails
            # phases is an int
            ret_phases = True
            phases = linspace(0., twopi, phases)
        rate, uncert = gl.rate_wb(self.lndfac, self.phases, self.phibins,
                       self.bins, self.offset, phases)
        # *** Is uncert a variance or std devn?
        if ret_phases:
            return phases, rate, uncert
        else:
            return rate, uncert
