"""Credibility interval HDR.

Author: Jose Vines

Taken from Pablo Cubillos MC3 and private communication with Monika Lendl
https://github.com/pcubillos/mc3/blob/master/mc3/stats/stats.py
"""

import numpy as np
from scipy.special import erf
from scipy.stats import gaussian_kde

def credibility_interval_hdr(dist, sigma=1.):
    """Calculate the highest density region for an empirical distribution.

    Reference: Hyndman, Rob J. 1996

    Parameters
    ----------
    dist: Array_like
        The posterior distribution for which the HDR is needed.
    sigma: float
        The confidence level in sigma notation. (e.g. 1 sigma = 68%)

    Returns
    -------
    best: float
        The value corresponding to the peak of the posterior distribution.
    low: float
        The minimum value of the HDR.
    high: float
        The maximum value of the HDR.

    Note: The HDR is capable of calculating more robust credible regions
    for multimodal distributions. It is identical to the usual probability
    regions of symmetric about the mean distributions. Using this then should
    lead to more realistic errorbars and 3-sigma intervals for multimodal
    distributions.

    """
    z = erf(sigma / np.sqrt(2))
    # First we estimate the PDF from the posterior distribution
    kde = gaussian_kde(dist)
    xmin, xmax = dist.min(), dist.max()
    xx = np.linspace(xmin, xmax, 1000)
    pdf = kde(xx)
    # This step is crucial else it doesn't work
    idx = np.argsort(pdf)[::-1]
    # Calculate the histogram
    hh, hx = np.histogram(dist, density=True, bins=1000)
    # Calculate the CDF
    cdf = np.cumsum(hh) * np.diff(hx)
    # Find where the CDF reaches 100*z%. Indices below this value are inside
    # The region corresponding to the confidence interval.
    idx_hdr = np.where(cdf >= z)[0][0]
    # Get the best value (i.e. the peak of the PDF)
    best = np.argmax(pdf)
    # Isolate the HDR
    hdr = pdf[idx][0:idx_hdr]
    # And get the minimum density
    hdr_min = hdr.min()
    # Min and max values are the extrema of the region and thus the confidence
    # interval
    low = np.min(xx[pdf > hdr_min])
    high = np.max(xx[pdf > hdr_min])
    return best, low, high