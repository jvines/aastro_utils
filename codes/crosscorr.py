import matplotlib.pyplot as plt
import scipy as sp
from PyAstronomy import pyasl
from scipy.optimize import curve_fit
from scipy.stats import norm


def ccf(wave, flux, mask='G2', rvmin=-300, rvmax=300, drv=0.1):
    """Compute the CCF of a spectrum.

    Based on PyAstronomy implementation of the cross-correlation.

    Parameters
    ----------
    wave: array_like
        The input wavelength array.
    flux: array_like
        The input flux array.
    mask: str
        The binary mask to use for the CCF. Options are G2, K0, K5 and M2
    rvmin: float, optional
        The minimum radial velocity to inspect in km/s. Default is -300 km/s
    rvmax: float, optional
        The maximum radial velocity to inspect in km/s. Default is 300 km/s
    drv: float, optional
        The radial velocity step. Default is 0.1 km/s

    Returns
    -------
    rvs: array_like
        The rv axis of the CCF.
    cc: array_like
        The median normalized CCF
    """
    # read mask, call crosscorr
    x1 = sp.arange(wave[0] - 200, wave[0], wave[1] - wave[0])
    x2 = sp.arange(wave[-1], wave[-1] + 200, wave[-1] - wave[-2])
    wtem = sp.hstack([x1, wave, x2])

    lines1, lines2, flux_l = sp.loadtxt('masks/' + mask + '.mas', unpack=True)
    ilines = sp.where((lines1 > wave[0]) & (lines2 < wave[-1]))[0]
    lines1_new = lines1[ilines]
    lines2_new = lines2[ilines]
    flux_l_new = flux_l[ilines]

    ftem = sp.zeros(wtem.size)

    for i in range(len(flux_l_new)):
        indices = sp.where((wtem >= lines1_new[i]) & (wtem <= lines2_new[i]))
        if indices[0].size > 0:
            ftem[indices[0]] = flux_l_new[i]
        del indices

    rvs, cc = pyasl.crosscorrRV(wave, flux, wtem, ftem, rvmin, rvmax, drv)
    return rvs, cc / sp.median(cc)


def ccf_plot(name, rvs, cc):
    """Plot of the CCF."""
    f, ax = plt.subplots(figsize=(12, 8))
    ax.plot(rvs, cc, lw=0.7, color='k')
    ax.set_ylabel('CCF')
    ax.set_xlabel('RV (km/s)')
    plt.savefig(name + '_ccf.pdf', bbox_inches='tight')


def ccf_fit(rvs, cc):
    """Quick gaussian fit of a CCF."""
    est = sp.argmin(cc)
    rv_min = rvs[est]
    cc_min = cc[est]
    p0 = [cc_min, rv_min, 1., 1.]
    popt, _ = curve_fit(gauss_fit, rvs, cc, p0=p0)
    return popt


def ccf_gaus_plot(name, rvs, cc, amp, rv, sig, off):
    """Plot of the CCF and the gaussian fit."""
    f, ax = plt.subplots(figsize=(12, 8))
    ax.scatter(rvs, cc, s=10, marker='.', color='k')
    ax.plot(rvs, gauss_fit(rvs, amp, rv, sig, off))
    ax.set_xlim([rv - 30, rv + 30])
    ax.set_ylabel('CCF')
    ax.set_xlabel('RV (km/s)')
    plt.savefig(name + '_ccf_fit.pdf', bbox_inches='tight')


def gauss_fit(r, a, mu, sig, off):
    """Gaussian model based on scipy.stats.norm"""
    return off + a * norm.pdf(r, loc=mu, scale=sig)
