from numpy import ones, zeros, shape, sum, newaxis, sqrt

# TODO: Use placeholder until vba is updated:
# import vba

__all__ = ['BA', 'BAObj']

# The placeholder:
vba = object()


# A global for alerting user to an untested method.
FIRSTTIME = True


class BA:
    """
    Implement the Bretthorst algorithm using a list of functions of the
    nonlinear parameters and an array of abscissas, a data array (abscissas and
    samples, with error sigmas as a 3rd dimension or as a separate, single sigma
    value), and an optional setup function.

    The setup function should take the nonlinear parameters as its arguments,
    and return an object that gets passed to each of the basis functions as
    its final argument.
    """

    def __init__(self, basis, data, sigma=None, setup=None):
        """Initialize:  Save basic info; calculate standardized data and d**2."""

#...    Basic member data.
        self.M = len(basis)
        self.N = len(data)
        self.basis = basis
        self.absc = data[:,0]
        self.smpls = data[:,1]
        if shape(data)[1] != 2 and sigma != None:
            raise ValueError("Data have sigmas but a sigma was provided!")
        if sigma != None:
            self.sigs = sigma*ones(self.N, float)
        elif (shape(data)[1] == 3):
            self.sigs = data[:,2]
        else:
            raise ValueError("Wrong dimensions for data!")
        self.setup = setup
        self.std_smpls = self.smpls / self.sigs
        self.dsqr = sum(self.std_smpls*self.std_smpls)

#...    We'll memoize the current params & results to avoid duplicating work.
        self.pars = None
        self.suf = None
        self.Q = None
        self.Jac = None

#...    Space for model values (aka the design matrix).
        self.modvals = zeros((self.M, self.N), float)

#...    We could also deal with constant terms here once for all, to handle
#...    purely linear terms.

    def marg_stats(self, *params):
        """
        Return the quadratic form, sufficient statistic, and Jacobian as a
        function of the passed nonlinear parameters.  These statistics allow
        calculation of the marginal density for the nonlinear parameters, and
        other inferences.

        Returns:
            Q - quadratic form (sum of squared residuals)
            suf - sufficient statistic
            Jac - Jacobian determinant

        Let dsqr be the sum of squares of the standardized sample values; then
        Q and suf are related, Q = dsqr - suf.

        The log marginal likelihood for the nonlinear parameters is
        -0.5*Q + log(J).
        """
        self.calc_marg_stats(*params)
        return self.Q, self.suf, self.Jac

    def calc_marg_stats(self, *params):
        """
        Calculate the marginal density and Jacobian vs. nonlinear parameters
        and memoize the results.
        """

        if (params == self.pars):
            return self.Q, self.suf, self.Jac

#...    Evaluate the standardized models on the data, first calling the setup
#...    function for the data if it is present.
        args = list(params)
        args.append(self.absc)
        args = tuple(args)
        if (self.setup != None):
            aux = self.setup(*args)
            args = list(args)
            args.append(aux)
            args = tuple(args)
        for a in range(self.M):
            g = self.basis[a]
            self.modvals[a,:] = g(*args) / self.sigs

#...    Calculate the metric and various derived quantities.  Store them
#...    for possible future use before returning.
        self.metric, self.L, self.Jac, self.proj, self.ampl, self.suf =\
            vba.metricAnalysis(self.modvals, self.std_smpls)
        self.Q = self.dsqr - self.suf
        self.pars = params
        self.Jac = 1. / self.Jac

    def amplitudes(self, *params):
        """Calculate the best-fit amplitudes."""

#...    If we haven't already done the metric calculations for these params,
#...    do them.
        if params != () and params != self.pars:
            # print "Calling margStats."
            self.calc_marg_stats(*params)

#...    Now just return the stored amplitudes..
        return self.ampl

    def covar(self):
        """Return the covariance matrix for the best-fit amplitudes."""

        global FIRSTTIME
        if FIRSTTIME:
            print('*** BA method covar is UNTESTED ***')  # Copied from old version--likely ok.
            FIRSTTIME = False

#...    If we haven't already done the metric calculations, complain!
        if self.pars == None:
            raise ValueError('Must calculate marginal first!')

        return vba.covar(self.L)

    def residuals(self, *params):
        """Calculate the standardized residuals."""

#...    If we haven't already done the metric calculations for these params,
#...    do them.
        if (params != () and params != self.pars):
            # print "Calling margStats."
            self.calc_marg_stats(*params)

#...    Subtract projection from data.
        model = sum(self.ampl[:,newaxis]*self.modvals, 0)
        return self.std_smpls - model


class BAObj:
    """
    Implement the Bretthorst algorithm using an object that provides basis
    function and data access through these methods and attributes:

        set_nonlin(*args) - set the nonlinear parameters
        std_smpls[i] - value of the i'th standardized sample
        std_basis[k] - method returning 1-D array of standardized predictions
            from basis function k when called as std_basis[k]()
    """

    def __init__(self, obj):
        """
        Initialize:  Specify number of basis functions, number of data samples,
        and the object providing basis and data access.
        """

#...    Basic member data.
        self.obj = obj
        self.M = len(obj.std_basis)
        self.N = len(obj.std_smpls)
        self.dsqr = sum(obj.std_smpls*obj.std_smpls)

#...    We'll memoize the current params & results to avoid duplicating work.
        self.pars = None
        self.suf = None
        self.Q = None
        self.Jac = None

#...    Space for model values (aka the design matrix).
        self.modvals = zeros((self.M, self.N), float)

    def marg_stats(self, *params):
        """
        Return the quadratic form, sufficient statistic, and Jacobian as a
        function of the passed nonlinear parameters.  These statistics allow
        calculation of the marginal density for the nonlinear parameters, and
        other inferences.

        Returns:
            Q - quadratic form (sum of squared residuals)
            suf - sufficient statistic
            Jac - Jacobian determinant

        Let dsqr be the sum of squares of the standardized sample values; then
        Q and suf are related, Q = dsqr - suf.

        If the prior for the amplitudes is taken to be flat and independent of
        the nonlinear parameters, the log marginal likelihood for the nonlinear
        parameters is
            -0.5*Q + log(Jac)
        However, if a subset of the basis functions becomes (nearly)
        degenerate for some choice of the nonlinear parameters, the likelihood
        becomes constant in one or more directions, and the marginal
        likelihood for that choice of nonlinear parameters will diverge
        because the amplitude integral diverges.  One way to guard against
        this is to make the amplitude prior *depend on the nonlinear
        parameters* proportionally to the volume element in sample space
        spanned by the basis.  Then nonlinear parameter values that cause
        "basis collapse" are penalized (the prior density vanishes for them).
        This corresponds to a prior for the amplitudes proportional to Jac, so
        the log marginal likelihood for the nonlinear parameters is then
            -0.5*Q
        """
        self.calc_marg_stats(*params)
        return self.Q, self.suf, self.Jac

    def calc_marg_stats(self, *params):
        """
        Calculate the marginal density and Jacobian vs. nonlinear parameters
        and memoize the results.
        """

        if (params == self.pars):
            return self.Q, self.suf, self.Jac

#...    Evaluate the standardized models on the data, first calling the setup
#...    function.
        self.obj.set_nonlin(*params)
        for a in range(self.M):
            self.modvals[a,:] = self.obj.std_basis[a]()

#...    Calculate the metric and various derived quantities.  Store them
#...    for possible future use before returning.
        self.metric, self.L, self.Jac, self.proj, self.ampl, self.suf =\
            vba.metricAnalysis(self.modvals, self.obj.std_smpls)
        self.Q = self.dsqr - self.suf
        self.pars = params
        self.Jac = 1. / self.Jac

    def amplitudes(self, *params):
        """Calculate the best-fit amplitudes."""

#...    If we haven't already done the metric calculations for these params,
#...    do them.
        if params != () and params != self.pars:
            # print "Calling margStats."
            self.calc_marg_stats(*params)

#...    Now just return the stored amplitudes.
        return self.ampl

    def covar(self):
        """Return the covariance matrix for the best-fit amplitudes."""

        global FIRSTTIME
        if FIRSTTIME:
            print('*** BA method covar is UNTESTED ***')  # Copied from old version--likely ok.
            FIRSTTIME = False

#...    If we haven't already done the metric calculations, complain!
        if self.pars == None:
            raise ValueError('Must calculate marginal first!')
        return vba.covar(self.L)

    def residuals(self, *params):
        """Calculate the standardized residuals."""

#...    If we haven't already done the metric calculations for these params,
#...    do them.
        if (params != () and params != self.pars):
            # print "Calling margStats."
            self.calc_marg_stats(*params)

#...    Subtract projection from data.
        model = sum(self.ampl[:,newaxis]*self.modvals, 0)
        return self.obj.std_smpls - model
