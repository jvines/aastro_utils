"""Cadence bin.

Author: Jose Vines
"""
import numpy as np


def cadence_bin(times, data, dt):
    """Bins timeseries data with cadence dt.

    Parameters:
    -----------
    times : array_like
        The times to bin in cadence dt.
    data : array_like
        Data corresponding to time times.
    dt : float
        Time cadence to bin into in minutes.

    Returns:
    --------
    binned_times : array_like
        The binned times
    binned_data : array_like
        The binned data corresponding to the median of all the original
        data values inside a bin.
    binned_errs : array_like
        The binned errors calculated as the square root of the variance of
        the data inside a given bin, divided by the square root of the
        number of data points inside the bin.

    """
    # First calculate the dt in JD
    dt *= 60 / 86400
    # Grab initial and final time in the timeseries
    ti = times[0]
    tf = times[-1]
    # Calculate number of bins
    n = int(np.ceil((tf - ti) / dt))
    binned_times = np.zeros(n - 1)
    binned_data = np.zeros(n - 1)
    binned_errs = np.zeros(n - 1)
    t = np.linspace(ti, tf, n)
    for i in range(0, n - 1):
        low = t[i] < times
        up = times < t[i + 1]
        bin_n = len(times[low * up])
        if ~np.any(low * up):
            continue
        binned_times[i] = np.median(times[low * up])
        binned_data[i] = np.median(data[low * up])
        binned_errs[i] = np.sqrt(np.var(data[low * up]) / bin_n)
    no_zeroes = binned_times != 0
    binned_times = binned_times[no_zeroes]
    binned_data = binned_data[no_zeroes]
    binned_errs = binned_errs[no_zeroes]
    return binned_times, binned_data, binned_errs
