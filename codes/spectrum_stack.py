"""Spectrum stack.

Stacks two or more echelle spectra for an object.

Author: Jose Vines
"""
import copy

import astropy.constants as const
import astropy.units as u
import scipy as sp
from astropy.io import fits
from tqdm import tqdm


def velocity_correction(instrument, spec_list):
    """Radial velocity correction for a list of spectra.

    Parameters
    ----------
    instrument : str
        The instrument with which the spectra was obtained.
    spec_list : array_like
        An array with filenames pointing to the fits files.

    """
    if instrument.lower() == 'feros':
        corrector = feros_velocity_correction
    for spec in tqdm(spec_list, desc='Spectra'):
        corrector(spec, True)
    pass


def feros_velocity_correction(spec, create_fits=False):
    """Radial velocity correction for a FEROS spectrum.

    Parameters
    ----------
    spec : file
        Fits file with FEROS spectra reduced with CERES.

    create_fits : bool, optional
        True to save a fits file with the result.

    Returns
    -------
    wavelength : array_like
        An array with the rest frame wavelengths.

    flux : array_like
        An array with the corresponding fluxes.

    """
    # Read fits file
    hdul = fits.open(spec)
    # Extract RV
    rv = hdul[0]['RV'] * u.km / u.s
    rv = rv.to(u.m / u. s)
    # Create gamma
    beta = rv / const.c
    gamma = 1 + beta.value
    # Extract wavelength per order
    wave = hdul[0].data[0, :, :]
    # Extract flux per order
    flux = hdul[0].data[5, :, 5]
    orders = wav.shape[0]
    wave_rest = copy.deepcopy(wave)
    # Move spectra to rest frame
    for o in orders:
        wave_rest[o, :] /= gamma
    if create_fits:
        # Create new fits file
        out = spec.split('.fits')[0] + '_rest_frame.fits'
        hdu = fits.PrimaryHDU(sp.stack(wave_rest, flux))
        hdu.writeto(out)
    hdul.close()
    return wave, flux


def median_combine(spec_list, nord, targ_name, ra, dec, plx):
    """Median combine rest frame spectra.

    Parameters
    ----------
    spec_list : array_like
        Array with spectra files.
    nord : int
        Number of echelle orders.
    targ_name : str
        Target's name.
    ra : float
        Target's right ascension in degrees.
    dec : float
        Target's declination in degrees.
    plx : float
        Target's parallax.

    """
    wavelengths = []
    fluxes = []
    for o in tqdm(range(nord), desc='Order'):
        combiner = []
        for spec in tqdm(spec_list, desc='Spectra'):
            hdul = fits.open(spec)
            flux = hdul[0].data[5, o, :]
            wave = hdul[0].data[0, o, :]
            combiner.append(flux)
            hdul.close()
        combiner = sp.array(combiner)
        combined = sp.median(combiner, axis=0)
        fluxes.append(combined)
        wavelengths.append(wave)
    final_waves = sp.vstack(wavelengths)
    final_fluxes = sp.vstack(fluxes)
    final_out = sp.stack([final_waves, final_fluxes])
    out = targ_name + 'median.fits'
    hdr = fits.Header()
    hdr['NAME'] = targ_name
    hdr['PLX'] = plx
    hdr['RA (deg)'] = ra
    hdr['DEC (deg)'] = dec
    hdu = fits.PrimaryHDU(final_out, header=hdr)
    hdu.writeto(out)
    pass
