"""get quadratures

Calculate the times of quadrature for a series of objects with given ephemerides
between two nights in a given observatory

Input file must be:

name ra(deg) dec(deg) epoch period

The output will be:
    A per target list containing the times of quadratures

"""

import argparse
from datetime import datetime
from datetime import timedelta

from astropy import units as u
from astropy.coordinates import AltAz
from astropy.coordinates import EarthLocation
from astropy.coordinates import SkyCoord
from astropy.coordinates import get_moon
from astropy.coordinates import get_sun
from astropy.time import Time
from astropy.utils.iers import conf as iers_conf

iers_conf.iers_auto_url = 'https://datacenter.iers.org/data/9/finals2000A.all'
mir = 'ftp://cddis.gsfc.nasa.gov/pub/products/iers/finals2000A.all'
iers_conf.iers_auto_url_mirror = mir


def arg_parse():
    """Parse command line arguments."""
    p = argparse.ArgumentParser()
    p.add_argument('infile', help='Path to file containing targets.')
    p.add_argument('n1', help='Night 1 in Y-m-d')
    p.add_argument('n2', help='Night 2 in Y-m-d')
    p.add_argument('observatory', help='Astropy name of the observatory')
    return p.parse_args()


def read_ephem_file(infile):
    """Read the ephem file."""
    name, ra, dec, epoch, period = [], [], [], [], []
    with open(infile, 'r') as f:
        for line in f:
            data = line.split()
            name.append(data[0])
            ra.append(float(data[1]))
            dec.append(float(data[2]))
            epoch.append(float(data[3]))
            period.append(float(data[4]))
    return name, ra, dec, epoch, period


def sun_is_down(time, observatory) -> bool:
    """Check if the Sun is below -14 deg altitude."""
    sun = get_sun(time).transform_to(AltAz(obstime=time, location=observatory))
    return sun.alt.value <= -14


def moon_is_away(time, ra, dec, observatory) -> bool:
    """Check if the moon is 30 deg away or more."""
    obj_coords = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs')
    moon = get_moon(time, location=observatory)
    sep = obj_coords.separation(moon).degree
    return sep >= 30


if __name__ == '__main__':
    args = arg_parse()
    # Observatory location.
    observatory = EarthLocation.of_site(args.observatory)
    # Read the ephem file.
    names, ras, decs, epochs, periods = read_ephem_file(args.infile)
    # Get times for the run.
    n1 = datetime.strptime(args.n1, '%Y-%m-%d') + timedelta(hours=12)
    # Add extra day so we have nights of each date.
    n2 = datetime.strptime(args.n2, '%Y-%m-%d') + timedelta(hours=36)
    n1_T = Time(n1, format='datetime', scale='utc', location=observatory)
    n2_T = Time(n2, format='datetime', scale='utc', location=observatory)
    # Loop over each object. Remember q1 is 0.25 phase and q2 is 0.75
    for name, ra, dec, epoch, period in zip(names, ras, decs, epochs, periods):
        print(f'Target: {name}')
        epoch_start = 0
        while epoch_start < n1_T.jd:
            epoch_start += period
        epoch_start -= period
        # Pull useable epochs
        current_epoch = epoch_start
        while current_epoch < n2_T.jd:
            q1 = current_epoch + 0.25 * period
            q2 = current_epoch - 0.25 * period
            q1_T = Time(q1, format='jd', scale='utc')
            q2_T = Time(q2, format='jd', scale='utc')
            if n1_T.jd <= q1 <= n2_T.jd:
                if sun_is_down(q1_T, observatory) and \
                        moon_is_away(q1_T, ra, dec, observatory):
                    print(f'\t{str(q1_T.datetime)[:16]}\t0.25')
            if n1_T.jd <= q2 <= n2_T.jd:
                if sun_is_down(q1_T, observatory) and \
                        moon_is_away(q1_T, ra, dec, observatory):
                    print(f'\t{str(q2_T.datetime)[:16]}\t0.75')
            current_epoch += period
