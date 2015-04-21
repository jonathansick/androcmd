#!/usr/bin/env python
# encoding: utf-8
"""
Pipelines for fitting PHAT CMDs.
"""
import os
from collections import namedtuple, OrderedDict

import numpy as np
from astropy.coordinates import Distance
import astropy.units as u
from astropy.table import Table

from m31hst import phat_v2_phot_path

from padova import AgeGridRequest
from padova.isocdata import join_isochrone_sets, Isochrone

from starfisher import ColorPlane
from starfisher import Lockfile
from starfisher.pipeline import (
    PipelineBase, IsochroneSetBase, DatasetBase, PlaneBase, LockBase)

PHAT_BANDS = ('F475W', 'F814W', 'F275W', 'F336W', 'F110W', 'F160W')
STARFISH = os.getenv("STARFISH")
Lim = namedtuple('Lim', 'x y')


class SolarZIsocs(IsochroneSetBase):
    """Solar metallicity isochrone set."""
    def __init__(self, **kwargs):
        self.isoc_args = kwargs.pop('isoc_args', dict())
        self.isoc_phases = kwargs.pop('isoc_phases', None)
        self.z_grid = [0.015, 0.019, 0.024]
        print "SolarZIsocs", kwargs
        super(SolarZIsocs, self).__init__(**kwargs)

    @property
    def bands(self):
        return ('F475W', 'F814W', 'F275W', 'F336W', 'F110W', 'F160W')

    @property
    def distance(self):
        return Distance(785. * u.kpc)

    def setup_isochrones(self):
        """Download Padova isochrones."""
        WFC3_BANDS = ['F275W1', 'F336W', 'F110W', 'F160W']
        ACS_BANDS = ['F475W', 'F814W']

        if self.isoc_args is None:
            self.isoc_args = {}
        if not os.path.exists(os.path.join(STARFISH, self.isoc_dir)):
            for z in self.z_grid:
                r_wfc3 = AgeGridRequest(z,
                                        min_log_age=6.6,
                                        max_log_age=10.13,
                                        delta_log_age=0.02,
                                        photsys='wfc3_wide', **self.isoc_args)
                r_acs = AgeGridRequest(z,
                                       min_log_age=6.6,
                                       max_log_age=10.13,
                                       delta_log_age=0.02,
                                       photsys='acs_wfc', **self.isoc_args)
                isoc_set = join_isochrone_sets(r_wfc3.isochrone_set,
                                               r_acs.isochrone_set,
                                               left_bands=WFC3_BANDS,
                                               right_bands=ACS_BANDS)
                for isoc in isoc_set:
                    isoc = Isochrone(isoc)
                    isoc.rename_column('F275W1', 'F275W')
                    if self.isoc_phases is not None:
                        sels = []
                        for p in self.isoc_phases:
                            sels.append(np.where(isoc['stage'] == p)[0])
                        s = np.concatenate(sels)
                        isoc = isoc[s]
                    isoc.export_for_starfish(os.path.join(STARFISH,
                                                          self.isoc_dir),
                                             bands=list(self.bands))
        super(SolarZIsocs, self).setup_isochrones()


class PhatCatalog(DatasetBase):
    """Mixin for PHAT photometry data.

    Photometry is lazy loaded to it is efficient to rebuild the pipeline
    object.
    """
    def __init__(self, **kwargs):
        self.brick = kwargs.pop('brick')
        self._phat_data = None
        super(PhatCatalog, self).__init__(**kwargs)

    def _load_phat_data(self):
        self._phat_data = Table.read(phat_v2_phot_path(self.brick),
                                     format='fits')
        phat_bands = ('F475W', 'F814W', 'F275W', 'F336W', 'F110W', 'F160W')
        # Normalize bandpass names
        for band in phat_bands:
            old_name = "_".join(band.lower, 'vega')
            self._phat_data.rename_column(old_name, band)

    def _make_phot_column(self, band):
        if not isinstance(band, basestring):
            band1, band2 = band
            return self._phat_data[band1] - self._phat_data[band2]
        else:
            return self._phat_data[band]

    def write_phot(self, x_band, y_band, data_root, suffix):
        if self._phat_data is None:
            self._load_phat_data()  # lazy loading

        x = self._make_phot_column(x_band)
        y = self._make_phot_column(y_band)
        phot_dtype = np.dtype([('x', np.float), ('y', np.float)])
        photdata = np.empty(len(self._phat_data), dtype=phot_dtype)
        photdata['x'][:] = x
        photdata['y'][:] = y

        path = data_root + suffix
        full_path = os.path.join(STARFISH, path)
        fit_dir = os.path.dirname(full_path)
        if not os.path.exists(fit_dir):
            os.makedirs(fit_dir)
        np.savetxt(full_path, photdata, delimiter=' ', fmt='%.4f')


class CompletePhatPlanes(PlaneBase):
    """Color planes for PHAT data."""
    def __init__(self, **kwargs):
        self._planes = OrderedDict([
            (('f475w', 'f814w'), make_f475w_f814w()),
            ('f475w_f814w_rgb', make_f475w_f814w_rgb()),
            (('f475w', 'f110w'), make_f475w_f110w()),
            (('f475w', 'f160w'), make_f475w_f160w()),
            (('f814w', 'f110w'), make_f814w_f110w()),
            (('f814w', 'f160w'), make_f814w_f160w()),
            (('f110w', 'f160w'), make_f110w_f160w()),
        ])
        print "CompletePhatPlanes", kwargs
        super(CompletePhatPlanes, self).__init__(**kwargs)

    @property
    def planes(self):
        return self._planes


class SolarLockfile(LockBase):
    """Lockfile mixin to create an iso-metallicity Hess set."""
    def __init__(self, **kwargs):
        print "SolarLockfile", kwargs
        super(SolarLockfile, self).__init__(**kwargs)

    def build_lockfile(self):
        self.lockfile = Lockfile(self.builder.read_isofile(), self.synth_dir,
                                 unbinned=False)

        # Bin young isochrones
        young_grid = np.linspace(6.5, 8.95, 10)
        for i, logage0 in enumerate(young_grid[:-1]):
            logage0 = logage0
            logage1 = young_grid[i + 1]
            z_str = "0019"
            mean_age = (logage0 + logage1) / 0.2
            name = "z{0}_{1:05.2f}".format(z_str, mean_age)
            self.lockfile.lock_box(name, (logage0, logage1), (0.014, 0.025))

        # Bin old isochrones
        old_grid = np.arange(1e9, 14 * 1e9, 1e9)
        for i, age0 in enumerate(old_grid[:-1]):
            logage0 = np.log10(age0 - 0.05 * 1e9)
            logage1 = np.log10(old_grid[i + 1])
            z_str = "0019"
            mean_age = (logage0 + logage1) / 0.2
            name = "z{0}_{1:05.2f}".format(z_str, mean_age)
            self.lockfile.lock_box(name, (logage0, logage1), (0.014, 0.025))


class SolarZPhatPipeline(CompletePhatPlanes, SolarZIsocs, PhatCatalog,
                         SolarLockfile, PipelineBase):
    """A pipeline for fitting PHAT bricks with solar metallicity isochrones."""
    def __init__(self, **kwargs):
        print "SolarZPhatPipeline", kwargs
        super(SolarZPhatPipeline, self).__init__(**kwargs)


def make_f475w_f814w(dpix=0.05, mag_lim=30.):
    lim = Lim(x=(-1, 5.), y=(25.5, 20.))
    plane = ColorPlane((PHAT_BANDS.index('F475W'),
                        PHAT_BANDS.index('F814W')),
                       PHAT_BANDS.index('F814W'),
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


def make_f475w_f814w_rgb(dpix=0.05, mag_lim=30.):
    lim = Lim(x=(1.2, 5.), y=(23.5, 20.))
    plane = ColorPlane((PHAT_BANDS.index('F475W'),
                        PHAT_BANDS.index('F814W')),
                       PHAT_BANDS.index('F814W'),
                       lim.x,
                       (min(lim.y), max(lim.y)),
                       mag_lim,
                       suffix='rgbopt',
                       x_label=r'$\mathrm{F475W}-\mathrm{F814W}$',
                       y_label=r'$\mathrm{F814W}$',
                       dpix=dpix,
                       nx=75)  # NOTE auto-calc failed to compute 75
    # plane.mask_region((3, 5), (28, 25))
    # plane.mask_region((3.5, 5), (25, 23))
    # plane.mask_region((4, 5), (23, 22.5))
    return plane


def make_f110w_f160w(dpix=0.05, mag_lim=30.):
    lim = Lim(x=(0.3, 1.3), y=(24., 16.5))
    plane = ColorPlane((PHAT_BANDS.index('F110W'),
                        PHAT_BANDS.index('F160W')),
                       PHAT_BANDS.index('F160W'),
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
    lim = Lim(x=(-0.8, 8.), y=(25., 17.5))
    plane = ColorPlane((PHAT_BANDS.index('F475W'),
                        PHAT_BANDS.index('F160W')),
                       PHAT_BANDS.index('F160W'),
                       lim.x,
                       (min(lim.y), max(lim.y)),
                       mag_lim,
                       suffix='f475f160',
                       x_label=r'$\mathrm{F475W}-\mathrm{F160W}$',
                       y_label=r'$\mathrm{F110W}$',
                       dpix=dpix)
    return plane


def make_f475w_f110w(dpix=0.05, mag_lim=30.):
    lim = Lim(x=(-0.8, 7.), y=(25., 18.))
    plane = ColorPlane((PHAT_BANDS.index('F475W'),
                        PHAT_BANDS.index('F110W')),
                       PHAT_BANDS.index('F110W'),
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
    plane = ColorPlane((PHAT_BANDS.index('F814W'),
                        PHAT_BANDS.index('F110W')),
                       PHAT_BANDS.index('F110W'),
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
    plane = ColorPlane((PHAT_BANDS.index('F814W'),
                        PHAT_BANDS.index('F160W')),
                       PHAT_BANDS.index('F160W'),
                       lim.x,
                       (min(lim.y), max(lim.y)),
                       mag_lim,
                       suffix='f814f160',
                       x_label=r'$\mathrm{F814W}-\mathrm{F160W}$',
                       y_label=r'$\mathrm{F160W}$',
                       dpix=dpix)
    return plane
