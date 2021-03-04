"""sample from distribution
Sample randomly from an empirical distribution.

Author: Jose Vines
"""
import numpy as np
from scipy.interpolate import interp1d


def sample_from_distribution(distribution, ini=None, end=None, size=100):
    """Take random samples from an empirical distribution.

    Parameters
    ----------
    distribution: array_like
        The empirical distribution from which we want to sample.
    ini: float, optional
        The starting point of the distribution.
    end: float, optional
        The ending point of the distribution.
    size: int, optional
        The number of samples we wish to extract.

    Returns
    -------
    samples: array_like
        The array with the random samples.
    """
    # First we calculate the CDF of the distribution
    h, hx = np.histogram(distribution, density=True, bins=499)
    cdf = np.zeros(500)  # ensure the first value of the CDF is 0
    cdf[1:] = np.cumsum(h) * np.diff(hx)
    # Now we interpolate the inverted cdf
    mn = ini if ini is not None else np.min(distribution)
    mx = end if end is not None else np.max(distribution)
    xx = np.linspace(mn, mx, 500)
    icdf = interp1d(cdf, xx)
    points = np.random.random_sample(size=size)
    return icdf(points)