#!/usr/bin/env python
# encoding: utf-8
"""
Pipeline for analysis of hte PHAT six-color photometry.
"""

import os
from glob import glob
from collections import namedtuple, OrderedDict

import numpy as np
from astropy.table import Table
from m31hst import phat_v2_phot_path

from astropy.coordinates import Distance
import astropy.units as u

from padova import AgeGridRequest
from padova.isocdata import join_isochrone_sets, Isochrone

from starfisher import ColorPlane
from starfisher import SimHess
from starfisher import LibraryBuilder
from starfisher import Lockfile
from starfisher import Synth
from starfisher import ExtinctionDistribution
from starfisher import ExtantCrowdingTable
from starfisher import SFH
from m31hst.phatast import PhatAstTable


STARFISH = os.getenv("STARFISH")
Lim = namedtuple('Lim', 'x y')
PHAT_BANDS = ['F275W', 'F336W', 'F475W', 'F814W', 'F110W', 'F160W']
WFC3_BANDS = ['F275W1', 'F336W', 'F110W', 'F160W']
ACS_BANDS = ['F475W', 'F814W']


class Pipeline(object):
    """Pipeline for Multi-CMD fitting and comparison"""
    def __init__(self, brick, root_dir, isoc_args=None, phases=None):
        super(Pipeline, self).__init__()
        self.brick = brick
        self.catalog = Catalog(brick)
        self.root_dir = root_dir
        self.isoc_dir = os.path.join(root_dir, 'isoc')
        self.lib_dir = os.path.join(root_dir, 'lib')
        self.synth_dir = os.path.join(root_dir, 'synth')

        self.z_grid = [0.015, 0.019, 0.024]

        self.get_isochrones(isoc_args=isoc_args, phases=phases)
        self.build_lockfile()
        self.planes = PhatPlanes()
        self.run_synth()

        self.fits = OrderedDict()

        self._solution_tables = {}

    def get_solution_table(self, key):
        if key not in self._solution_tables:
            tbl = self.fits[key].solution_table()
            self._solution_tables[key] = tbl
        return tbl

    def get_isochrones(self, isoc_args=None, phases=None):
        if isoc_args is None:
            isoc_args = {}
        if not os.path.exists(os.path.join(STARFISH, self.isoc_dir)):
            for z in self.z_grid:
                r_wfc3 = AgeGridRequest(z,
                                        min_log_age=6.6,
                                        max_log_age=10.13,
                                        delta_log_age=0.02,
                                        photsys='wfc3_wide', **isoc_args)
                r_acs = AgeGridRequest(z,
                                       min_log_age=6.6,
                                       max_log_age=10.13,
                                       delta_log_age=0.02,
                                       photsys='acs_wfc', **isoc_args)
                isoc_set = join_isochrone_sets(r_wfc3.isochrone_set,
                                               r_acs.isochrone_set,
                                               left_bands=WFC3_BANDS,
                                               right_bands=ACS_BANDS)
                for isoc in isoc_set:
                    isoc = Isochrone(isoc)
                    isoc.rename_column('F275W1', 'F275W')
                    if phases is not None:
                        sels = []
                        for p in phases:
                            sels.append(np.where(isoc['stage'] == p)[0])
                        s = np.concatenate(sels)
                        isoc = isoc[s]
                    isoc.export_for_starfish(os.path.join(STARFISH,
                                                          self.isoc_dir),
                                             bands=PHAT_BANDS)

        d = Distance(785 * u.kpc)
        self.builder = LibraryBuilder(self.isoc_dir, self.lib_dir,
                                      nmag=len(PHAT_BANDS),
                                      dmod=d.distmod.value,
                                      iverb=3)
        if not os.path.exists(self.builder.full_isofile_path):
            self.builder.install()

    def build_lockfile(self):
        if not os.path.exists(os.path.join(STARFISH, self.synth_dir)):
            os.makedirs(os.path.join(STARFISH, self.synth_dir))

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

    def run_synth(self):
        full_synth_dir = os.path.join(STARFISH, self.synth_dir)
        if not os.path.exists(full_synth_dir):
            os.makedirs(full_synth_dir)

        # Use PHAT AST from the outer field (field 0)
        crowd_path = os.path.join(self.synth_dir, "crowding.dat")
        full_crowd_path = os.path.join(STARFISH, crowd_path)
        tbl = PhatAstTable()
        tbl.write_crowdfile_for_field(full_crowd_path, 0,
                                      bands=PHAT_BANDS)
        crowd = ExtantCrowdingTable(crowd_path)

        # No extinction, yet
        young_av = ExtinctionDistribution()
        old_av = ExtinctionDistribution()
        rel_extinction = np.ones(len(PHAT_BANDS), dtype=float)
        for av in (young_av, old_av):
            av.set_uniform(0.)

        self.synth = Synth(self.synth_dir, self.builder, self.lockfile, crowd,
                           rel_extinction,
                           young_extinction=young_av,
                           old_extinction=old_av,
                           planes=self.planes.all_planes,
                           mass_span=(0.08, 150.),
                           nstars=10000000)
        if len(glob(os.path.join(STARFISH, self.synth_dir, "z*"))) == 0:
            self.synth.run_synth(n_cpu=4)

    def fit_planes(self, key, color_planes, phot_colors):
        fit_dir = os.path.join(self.root_dir, key)
        data_root = os.path.join(fit_dir, "phot.")
        for band1, band2 in phot_colors:
            print data_root, band1, band2
            self.catalog.write(band1, band2, data_root)
        sfh = SFH(data_root, self.synth, fit_dir, planes=color_planes)
        if not os.path.exists(sfh.full_outfile_path):
            sfh.run_sfh()
        self.fits[key] = sfh


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
        full_path = os.path.join(STARFISH, path)
        fit_dir = os.path.dirname(full_path)
        if not os.path.exists(fit_dir):
            os.makedirs(fit_dir)
        np.savetxt(full_path, photdata, delimiter=' ', fmt='%.4f')


class PhatPlanes(object):
    """Color planes for PHAT data."""
    def __init__(self):
        super(PhatPlanes, self).__init__()
        self._planes = OrderedDict([
            (('f475w', 'f814w'), make_f475w_f814w()),
            (('f475w', 'f110w'), make_f475w_f110w()),
            (('f475w', 'f160w'), make_f475w_f160w()),
            (('f814w', 'f110w'), make_f814w_f110w()),
            (('f814w', 'f160w'), make_f814w_f160w()),
            (('f110w', 'f160w'), make_f110w_f160w()),
        ])

        self._sim_hess_planes = {}

    def get(self, band1, band2):
        return self._planes[(band1, band2)]

    def get_sim_hess(self, band1, band2, synth, lockfile):
        key = (band1, band2)
        if key not in self._sim_hess_planes:
            sh = SimHess(synth, self._planes[key],
                         np.ones(len(lockfile.active_groups)))
            self._sim_hess_planes[key] = sh
        return self._sim_hess_planes[key]

    @property
    def all_planes(self):
        return [p for k, p in self._planes.iteritems()]


def make_f475w_f814w(dpix=0.05, mag_lim=30.):
    lim = Lim(x=(-1, 5), y=(25.5, 20))
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


def make_f110w_f160w(dpix=0.05, mag_lim=30.):
    lim = Lim(x=(0.3, 1.3), y=(24, 16.5))
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
    lim = Lim(x=(-0.8, 8), y=(25, 17.5))
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
    lim = Lim(x=(-0.8, 7), y=(25, 18))
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
