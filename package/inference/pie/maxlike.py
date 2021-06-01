from .inference_base import Inference, maximize


class MaxLikeInference(Inference):
    """
    Base class for maximum likelihood inferences.
    """

    # Larger values of likelihood are preferred!
    extremize = maximize

    # Set default minimization method and tolerance
    min_method = 'powell'
    min_tol = 1.e-6

    def score(self):
        """
        The function to optimize in a fit---here the log likelihood.
        """
        self._log_like = 0.
        for name in self.predictor_names:
            self._log_like += getattr(self, name).log_like()
        return self._log_like
