"""Credibility interval.

Author: Jose Vines
"""
import scipy as sp


def credibility_interval(post, alpha=1.):
    """Calculate bayesian credibility interval.

    Parameters:
    -----------
    post : array_like
        The posterior sample over which to calculate the bayesian credibility
        interval.
    alpha : float, optional
        Confidence level.
    Returns:
    --------
    med : float
        Median of the posterior.
    low : float
        Lower part of the credibility interval.
    up : float
        Upper part of the credibility interval.

    """
    z = erf(alpha / sp.sqrt(2))

    lower_percentile = 100 * (1 - z) / 2
    upper_percentile = 100 * (1 + z) / 2
    low, med, up = sp.percentile(
        post, [lower_percentile, 50, upper_percentile]
    )
    return med, low, up
