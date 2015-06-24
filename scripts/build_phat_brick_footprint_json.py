#!/usr/bin/env python
# encoding: utf-8
"""
Build a json file of PHAT footprints for plotting.

The footprints are available from XXX.

2015-06-24 - Created by Jonathan Sick
"""

import json

import astropy.io.fits as fits
from astropy.wcs import WCS
from m31hst import phat_brick_path


def main():
    footprints = {}
    bricks = range(1, 24)
    _ = zip(bricks, map(phat_brick_path, bricks, ['f814w'] * len(bricks)))
    for b, path in _:
        with fits.open(path) as f:
            header = f[1].header
            wcs = WCS(header)
            footprint = wcs.calc_footprint()
            footprints[b] = footprint.tolist()
    with open('phat_brick_footprints.json', 'w') as f:
        json.dump(footprints, f, indent=4)


if __name__ == '__main__':
    main()
