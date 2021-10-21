"""Credibility interval.

Author: Jose Vines
"""
import numpy as np
from scipy.special import erf


def credibility_interval(post, alpha=1., axis=None):
    """Calculate bayesian credibility interval.

    Parameters:
    -----------
    post: array_like
        The posterior sample over which to calculate the bayesian credibility
        interval.
    alpha: float, optional
        Confidence level.
    axis: int, optional
        Along which axis will the CI be calculated.

    Returns:
    --------
    med: float
        Median of the posterior.
    low: float
        Lower part of the credibility interval.
    up: float
        Upper part of the credibility interval.

    """
    z = erf(alpha / np.sqrt(2))

    lower_percentile = 100 * (1 - z) / 2
    upper_percentile = 100 * (1 + z) / 2
    low, med, up = np.percentile(
        post, [lower_percentile, 50, upper_percentile], axis=axis
    )
    return med, low, up
