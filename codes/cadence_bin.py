"""
Cadence bin.

Author: Jose Vines
"""


def cadence_bin(times, data, dt):
    """Bins timeseries data with cadence dt.

    Parameters:
    -----------
    times : array_like
        The times to bin in cadence dt
    data : array_like
        Data corresponding to time times
    dt : float
        Time cadence to bin into

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
    def cadence_bin(times, data, dt):
        i = 0
        ti = times[0]
        tf = times[-1]
        n = int(sp.ceil((tf - ti) / dt))
        binned_times = sp.zeros(n - 1)
        binned_data = sp.zeros(n - 1)
        binned_errs = sp.zeros(n - 1)
        t = sp.linspace(ti, tf, n)
        for i in range(0, n - 1):
            low = t[i] < times
            up = times < t[i + 1]
            bin_n = len(times[low * up])
            if ~sp.any(low * up):
                continue
            binned_times[i] = sp.median(times[low * up])
            binned_data[i] = sp.median(data[low * up])
            binned_errs[i] = sp.sqrt(sp.var(data[low * up]) / bin_n)
        return binned_times, binned_data, binned_errs
