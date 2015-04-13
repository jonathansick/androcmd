#!/usr/bin/env python
# encoding: utf-8
"""
Pipeline for analysis of hte PHAT six-color photometry.
"""

import os
from collections import namedtuple

import numpy as np
from astropy.table import Table
from m31hst import phat_v2_phot_path

from starfisher import ColorPlane


Lim = namedtuple('Lim', 'x y')
wfc3_bands = ['F475W', 'F814W', 'F110W', 'F160W']


class Catalog(object):
    """Brick data catalog."""
    def __init__(self, brick):
        super(Catalog, self).__init__()
        self.data = Table.read(phat_v2_phot_path(brick), format='fits')

    def written_path(self, band1, band2, data_root):
        """Path to the output data file."""
        ext = '{0}{1}'.format(band1.rstrip('w'), band2.rstrip('w'))
        path = data_root + ext
        return path

    def write(self, band1, band2, data_root):
        """Write a band1-band2 vs band2 photometry catalog."""
        bands = (band1, band2)
        keys = ['{0}_vega'.format(band) for band in bands]
        phot_dtype = np.dtype([('x', np.float), ('y', np.float)])
        photdata = np.empty(len(self.data), dtype=phot_dtype)
        photdata['x'][:] = self.data[keys[0]] - self.data[keys[1]]
        photdata['y'][:] = self.data[keys[1]]

        path = self.written_path(band1, band2, data_root)
        fit_dir = os.path.dirname(path)
        if not os.path.exists(fit_dir):
            os.makedirs(fit_dir)
        np.savetxt(path, photdata, delimiter=' ', fmt='%.4f')


class PhatPlanes(object):
    """Color planes for PHAT data."""
    def __init__(self):
        super(PhatPlanes, self).__init__()
        self._planes = {
            ('f475w', 'f814w'): make_f475w_f814w(),
        }


def make_f475w_f814w(dpix=0.05, mag_lim=30.):
    lim = Lim(x=(-1, 5), y=(25.5, 20))
    plane = ColorPlane((wfc3_bands.index('F475W'),
                        wfc3_bands.index('F814W')),
                       wfc3_bands.index('F814W'),
                       lim.x,
                       (min(lim.y), max(lim.y)),
                       mag_lim,
                       suffix='f475f814',
                       x_label=r'$\mathrm{F475W}-\mathrm{F814W}$',
                       y_label=r'$\mathrm{F814W}$',
                       dpix=dpix)
    plane.mask_region((3, 5), (28, 25))
    plane.mask_region((3.5, 5), (25, 23))
    plane.mask_region((4, 5), (23, 22.5))
    return plane


def make_f110w_f160w(dpix=0.05, mag_lim=30.):
    lim = Lim(x=(0.3, 1.3), y=(24, 16.5))
    plane = ColorPlane((wfc3_bands.index('F110W'),
                        wfc3_bands.index('F160W')),
                       wfc3_bands.index('F160W'),
                       lim.x,
                       (min(lim.y), max(lim.y)),
                       mag_lim,
                       suffix='f110f160',
                       x_label=r'$\mathrm{F110W}-\mathrm{F160W}$',
                       y_label=r'$\mathrm{F160W}$',
                       dpix=dpix)
    plane.mask_region((-1., 0.), (22., 16))
    plane.mask_region((0, 0.3), (22., 16))
    plane.mask_region((0.3, 0.7), (20., 16))
    plane.mask_region((0.7, 0.8), (19., 16))
    plane.mask_region((0.8, 0.9), (18., 16))
    plane.mask_region((1.1, 1.5), (28, 21))
    return plane


def make_f475w_f160w(dpix=0.05, mag_lim=30.):
    lim = Lim(x=(-0.8, 8), y=(25, 17.5))
    plane = ColorPlane((wfc3_bands.index('F475W'),
                        wfc3_bands.index('F160W')),
                       wfc3_bands.index('F160W'),
                       lim.x,
                       (min(lim.y), max(lim.y)),
                       mag_lim,
                       suffix='f475f160',
                       x_label=r'$\mathrm{F475W}-\mathrm{F160W}$',
                       y_label=r'$\mathrm{F110W}$',
                       dpix=dpix)
    return plane


def make_f475w_f110w(dpix=0.05, mag_lim=30.):
    lim = Lim(x=(-0.8, 7), y=(25, 18))
    plane = ColorPlane((wfc3_bands.index('F475W'),
                        wfc3_bands.index('F110W')),
                       wfc3_bands.index('F110W'),
                       lim.x,
                       (min(lim.y), max(lim.y)),
                       mag_lim,
                       suffix='f475f110',
                       x_label=r'$\mathrm{F475W}-\mathrm{F110W}$',
                       y_label=r'$\mathrm{F110W}$',
                       dpix=dpix)
    return plane


def make_f814w_f110w(dpix=0.05, mag_lim=30.):
    lim = Lim(x=(-0.1, 1.8), y=(25, 19))
    plane = ColorPlane((wfc3_bands.index('F814W'),
                        wfc3_bands.index('F110W')),
                       wfc3_bands.index('F110W'),
                       lim.x,
                       (min(lim.y), max(lim.y)),
                       mag_lim,
                       suffix='f814f110',
                       x_label=r'$\mathrm{F814W}-\mathrm{F110W}$',
                       y_label=r'$\mathrm{F110W}$',
                       dpix=dpix)
    return plane


def make_f814w_f160w(dpix=0.05, mag_lim=30.):
    lim = Lim(x=(-0.5, 3), y=(24, 17.5))
    plane = ColorPlane((wfc3_bands.index('F814W'),
                        wfc3_bands.index('F160W')),
                       wfc3_bands.index('F160W'),
                       lim.x,
                       (min(lim.y), max(lim.y)),
                       mag_lim,
                       suffix='f814f160',
                       x_label=r'$\mathrm{F814W}-\mathrm{F160W}$',
                       y_label=r'$\mathrm{F160W}$',
                       dpix=dpix)
    return plane
