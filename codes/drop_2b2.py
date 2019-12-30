#!/usr/bin/env python
"""Remove spectra with 2x2 binning.

Author: Jose Vines
"""
import argparse
import glob
import os

from astropy.io import fits
from tqdm import tqdm

if __name__ == '__main__':
    description = '''
    Removes FEROS RAW spectra with 2x2 binning. Input is the directory
    holding the raw spectra.
    '''
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('in_dir', metavar='d', type=str,
                        help='FEROS spectra directory.')

    args = parser.parse_args()
    dir_to_wave = args.in_dir

    files = glob.glob(dir_to_wave + '*.fits')

    print('Reading spectra')

    for f in tqdm(files):
        h = fits.open(f)[0].header
        if h['CDELT1'] == 2. or h['CDELT2'] == 2:
            print('Removed ' + f)
            os.remove(f)
