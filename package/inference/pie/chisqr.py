from .inference_base import Inference, minimize


class ChisqrInference(Inference):
    """Base class for inferences made by minimizing a chi**2 statistic
    (a weighted sum of squared residuals)."""

    # Smaller values of chi**2 are preferred!
    extremize = minimize

    # Set default minimization method and tolerance
    min_method = 'powell'
    min_tol = 1.e-6

    def score(self):
        """The function to *minimize* in a fit; i.e., parameters with
        *smaller* values of this function are preferred to those with
        larger values.

        For chi**2 fitting, this is simply chi**2."""
        chi2 = 0.
        for name in self.predictor_names:
            chi2 += getattr(self, name).chisqr()
        return chi2
