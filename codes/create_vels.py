#!/Users/jvines/anaconda3/bin/python
"""
create_vels.py

Author: Jose Vines
Creates or update a .dat file, adding a new velocity, bisector, time, etc.
If the velocity is already present, it doesn't add it.
"""

from __future__ import print_function

import argparse
import os

import numpy as np


def load_vels(fname, unpack=False):
    return np.loadtxt(fname, unpack=unpack)


if __name__ == '__main__':
    description = '''
    Creates or updates vels files from the results of CERES pipeline
    '''
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('read_from_file', metavar='f', type=str,
                        help='CERES results.')
    parser.add_argument('vels_dir', metavar='v', type=str,
                        help='vels directory.')

    args = parser.parse_args()
    dir_to_res = args.read_from_file + 'results.txt'
    dir_to_vels = args.vels_dir

    # READ OR CREATE VELS FILE!!!

    print('Reading results file')
    data = np.loadtxt(dir_to_res, dtype=object)
    n_vels = len(data.shape)

    for i in range(data.shape[0]):
        skip = False
        if n_vels > 1:
            name = data[i, 0]
            info = data[i, 1:6]  # JD, RV, RV_e, BIS, BIS_e
            instr = data[i, 6]
        else:
            name = data[0]
            info = data[1:6]
            instr = data[6]
        info = info.astype(float)
        fname = dir_to_vels + name + '_' + instr.upper() + '.dat'
        has_content = False
        exists = os.path.isfile(fname)
        if exists:
            has_content = os.stat(fname).st_size

        if exists and has_content:
            print('Updating ' + fname)
            vels = load_vels(fname)
            jd = float(info[0])
            if jd not in vels:
                new_vels = np.vstack([info, vels])
                idx = np.argsort(new_vels, axis=0)[:, 0]
                np.savetxt(fname, new_vels[idx], delimiter=' ', fmt='%s')
            else:
                print('Velocity data already in file. Skipping.')
                skip = True
        elif not skip:
            print('Writing ' + fname)
            np.savetxt(fname, info.reshape(1, info.shape[0]),
                       delimiter=' ', fmt='%s')
        if n_vels == 1:
            break
