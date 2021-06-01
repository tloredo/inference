from .inference_base import Inference, maximize
from ..deriv import dllsteps, obsinfo

# Conceptually, it's appropriate to subclass this from MaxLikeInference, but
# there is so little inherited that it doesn't make any sense to subclass.


class BayesianInference(Inference):
    """Base class for Bayesian inferences made using the product of
    a prior distribution for parameters and a likelihood function."""

    # Larger values of prior*likelihood are preferred!
    extremize = maximize

    # Set default minimization method and tolerance
    min_method = 'powell'
    min_tol = 1.e-6

    # There is intentionally no default MCMC method; a good method
    # must be problem-specific.
    MCMC_method = None

    # def prior(self):
    # We could put a prior/logprior prototype here for users to override,
    # but if they define a prior in their SignalModel classes, it won't
    # override the prototype if BayesianInference is the first base class
    # in their inference class.

    def score(self):
        """
        The function to optimize in a fit---here the log(prior*likelihood).

        Note that a rigorous Bayesian 'optimum' requires specification of
        a utility (or loss) function; the decision-theoretic optimum maximizes
        the expected utility.
        """
        self._log_like = 0.
        for name in self.predictor_names:
            self._log_like += getattr(self, name).log_like()
        self._log_prior = self.log_prior()
        return self._log_like + self._log_prior

    def logp(self):
        """
        The log 'posterior' (i.e., prior*likelihood, unnormalized) for the
        current parameter values.
        """
        return self.score()

    def _logp_func(self, args):
        """
        The log posterior (prior*likelihood), implemented as a function of
        an array (or other sequence) of values of varying parameters.

        This is intended for use by algorithms that require a function with a
        sequence argument, e.g., the observed information algorithm.
        """
        # *** Note this changes *all* varying parameters, even those whose
        # actual value is not changed.  See the note for _score in inference.py.
        for param, val in zip(self.varying_params, args):
            param.set_value(val)
        return self.score()

    def obs_info(self, steps=None, dlogp=1., n=3, err=False):
        """
        Calculated the observed information matrix for varying parameters,
        assuming their current values correspond to a mode (at which a fit
        has already been evaluated).
        """
        # *** Expand docstring.
        # Treat the current state as the mode.
        logp_max = self._log_like + self._log_prior
        mode = []
        for param in self.varying_params:
            mode.append(param.get_value())
        # Find scales for the parameters that change logp by dlogp.
        # Start with the steps arg if provided; otherwise use current
        # varying parameter scales.
        if steps is None:
            deltas = []
            for param in self.varying_params:
                deltas.append(param.delta)
        else:
            deltas = steps
        deltas = dllsteps(self._logp_func, mode, deltas, dlogp, logp_max)
        # Use those step sizes for estimating the observed info matrix.
        info, ierr = obsinfo(self._logp_func, mode, deltas, logp_max, niter=n)
        if err:
            return info, ierr
        else:
            return info

    def marg(self):
        """
        Return the log integral of prior*likelihood, with the integral
        performed over varying parameters via adaptive quadrature.
        """
        raise NotImplementedError

    def laplace(self):
        """
        Return an approximation to the log of the integral of the
        prior*likelihood over varying parameters, calculated using the Laplace 
        approximation, i.e., by maximizing and multiplying by the determinant
        of the information matrix.
        """
        # *** Enable this on a grid.
        raise NotImplementedError

    def set_proposer(self):
        """Define the proposal distribution for MCMC posterior sampling."""
        raise NotImplementedError

    def proposal(self):
        """Return a candidate new state for the Markov chain."""
        raise NotImplementedError

    def MH_step(self):
        """Perform a single step of the Metropolis-Hastings sampler."""
        raise NotImplementedError

    def MH_steps(self, n, thin=None):
        """Perform n steps of the Metropolis-Hastings sampler.  If thin=t is
        provided, store a thinned version of the resulting time series,
        keeping only every t'th sample."""
        raise NotImplementedError
