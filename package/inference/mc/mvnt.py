from numpy import zeros, ones, diag, exp, random
from scipy.stats import randint
import _mvnt


class MVN(object):
    """Multivariate normal distribution:  Sampler and PDF."""

    def __init__(self, mu, covar=None, icovar=None):
        self.mu = mu
        self.ndim = len(mu)
        if covar is None and icovar is None:  # unit covar if unspecified
            covar = diag(ones(self.ndim))
        if covar != None:
            if icovar != None:
                raise ValueError('Specify only one of covar/icovar!')
            self.covar = covar
            ccopy, self.L, self.norm, err = _mvnt.mvnini(covar, 0)
        else:
            self.icovar = icovar
            self.covar, self.L, self.norm, err = _mvnt.mvnini(covar, 1)
        if err != 0: raise RuntimeError('Linear algebra failure in MVN init!')
        self.work = zeros(self.ndim, float)  # workspace for density calc'n

    def sample(self, n=None):
        """Return a multivariate normal sample."""
        if n is None:
            snsamps = random.standard_normal(self.ndim)
            return _mvnt.mvnsamp(self.mu, self.L, snsamps)
        samps = zeros((n, self.ndim), float)
        for i in range(n):
            snsamps = random.standard_normal(self.ndim)
            samps[i] = _mvnt.mvnsamp(self.mu, self.L, snsamps)
        return samps

    def sample_q(self, n=None):
        """Return a multivariate normal sample and the value of its
        associated quadratic form."""
        if n is None:
            snsamps = random.standard_normal(self.ndim)
            return _mvnt.mvnsampq(self.mu, self.L, snsamps)
        samps = zeros((n, self.ndim), float)
        qvals = zeros(n, float)
        for i in range(n):
            snsamps = random.standard_normal(self.ndim)
            samps[i], qvals[i] = _mvnt.mvnsampq(self.mu, self.L, snsamps)
        return samps, qvals

    def pdf(self, samp):
        """Multivariate normal probability density for samp."""
        return _mvnt.mvnden(samp, self.mu, self.norm, self.L, self.work)

    def pdfq(self, Q):
        """Multivariate normal probability density for a sample with
        quadratic form equal to Q."""
        return self.norm*exp(-Q/2.)


class MVT(object):
    """Multivariate Student's t distribution:  Sampler and PDF."""

    def __init__(self, nu, mu, hess):
        self.nu = nu
        self.hnu = 0.5*nu
        self.mu = mu
        self.ndim = len(mu)
        self.hess = hess
        self.ihess, self.L, self.norm, err = _mvnt.mvtini(nu, hess)
        if err != 0: raise RuntimeError('Linear algebra failure in MVT init!')
        self.work = zeros(self.ndim, float)  # workspace for density calc'n

    def sample(self, n=None):
        """Return a multivariate t sample."""
        if n is None:
            snsamps = random.standard_normal(self.ndim)
            gamsamp = random.gamma(self.hnu)
            return _mvnt.mvtsamp(self.mu, self.L, self.nu, snsamps, gamsamp)
        samps = zeros((n, self.ndim), float)
        for i in range(n):
            snsamps = random.standard_normal(self.ndim)
            gamsamp = random.gamma(self.hnu)
            samps[i] = _mvnt.mvtsamp(self.mu, self.L, self.nu, snsamps, gamsamp)
        return samps

    def sample_q(self):
        """Return a multivariate t sample and the value of its
        associated quadratic form."""
        if n is None:
            snsamps = random.standard_normal(self.ndim)
            gamsamp = random.gamma(self.hnu)
            return _mvnt.mvtsampq(self.mu, self.L, self.nu, snsamps, gamsamp)
        samps = zeros((n, self.ndim), float)
        qvals = zeros(n, float)
        for i in range(n):
            snsamps = random.standard_normal(self.ndim)
            gamsamp = random.gamma(self.hnu)
            samps[i], qvals[i] = _mvnt.mvtsampq(self.mu, self.L, self.nu, snsamps, gamsamp)
        return samps

    def pdf(self, samp):
        """Multivariate t probability density for samp."""
        return _mvnt.mvtden(samp, self.mu, self.nu, self.norm, self.L, self.work)

    def pdfq(self, Q):
        """Multivariate t probability density for a sample with quadratic
        form equal to Q."""
        return _mvnt.mvtdenq(self.ndim, self.nu, self.norm, Q)


class MVNKDE(object):
    """Multivariate normal kernel density estimate:  Sampler and PDF."""

    def __init__(self, nodes, scale, covar=None, icovar=None):
        if covar is None and icovar is None:
            raise ValueError('Must specify covar or icovar!')
        self.nodes = nodes.copy()
        self.tnodes = nodes.transpose().copy()
        self.nnodes, self.ndim = nodes.shape
        self.scale = scale
        if covar != None:
            if icovar != None:
                raise ValueError('Specify only one of covar/icovar!')
            self.covar = covar
            ccopy, self.L, self.norm, err = _mvnt.mvnini(covar, 0)
        else:
            self.icovar = icovar
            self.covar, self.L, self.norm, err = _mvnt.mvnini(covar, 1)
        if err != 0: raise RuntimeError('Linear algebra failure in MVN init!')
        self.work = zeros(self.ndim, float)  # workspace for density calc'n
        self.origin = zeros(self.ndim, float)
        self.randint = randint(0, self.nnodes)

    def sample(self, n=None):
        """
        Return one or more samples from the KDE.
        Randomly choose a component, then return a scaled MVN sample from
        that component.
        """
        if n == None:
            n = self.randint.rvs()
            snsamps = random.standard_normal(self.ndim)
            displ = _mvnt.mvnsamp(self.origin, self.L, snsamps)
            return self.nodes[n] + self.scale*displ
        samps = zeros((n, self.ndim), float)
        for i in range(n):
            n = self.randint.rvs()
            snsamps = random.standard_normal(self.ndim)
            displ = _mvnt.mvnsamp(self.origin, self.L, snsamps)
            samps[i] = self.nodes[n] + self.scale*displ
        return samps

    def pdf(self, samp):
        """Multivariate normal probability density for samp."""
        return _mvnt.mvnkde(samp, self.tnodes, self.norm, self.L, self.work,
                            self.scale)

    def set_node(self, n, pt):
        """
        Change the location of node n.
        """
        self.nodes[n] = pt[:]
        self.tnodes[:,n] = pt[:]

