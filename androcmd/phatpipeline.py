#!/usr/bin/env python
# encoding: utf-8
"""
Pipelines for fitting PHAT CMDs.
"""
import os
from glob import glob
from functools import partial

import matplotlib as mpl
import matplotlib.pyplot as plt
import palettable

import numpy as np
import astropy
from astropy.coordinates import Distance
import astropy.units as u
from astropy.table import Table

from m31hst import phat_v2_phot_path, phat_brick_path
from m31hst.phatast import PhatAstTable

from padova import AgeGridRequest
from padova.isocdata import join_isochrone_sets, Isochrone

from starfisher import Lockfile
from starfisher import ExtantCrowdingTable
from starfisher.dust import ExtinctionDistribution
from starfisher.pipeline import (
    PipelineBase, IsochroneSetBase, DatasetBase, LockBase,
    CrowdingBase, ExtinctionBase)

from androcmd.planes import RgbPhatPlanes, CompletePhatPlanes
from androcmd.dust import mw_Av, phat_rel_extinction, LewisDustLaw


PHAT_BANDS = ('F475W', 'F814W', 'F275W', 'F336W', 'F110W', 'F160W')
STARFISH = os.getenv("STARFISH")


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
        print "Running SolarZIsocs setup_isochrones"
        WFC3_BANDS = ['F275W1', 'F336W', 'F110W', 'F160W']
        ACS_BANDS = ['F475W', 'F814W']

        if not os.path.exists(os.path.join(STARFISH, self.isoc_dir)):
            make_isocs = True
        elif len(glob(os.path.join(STARFISH, self.isoc_dir, "*"))) == 0:
            make_isocs = True
        else:
            make_isocs = False

        if self.isoc_args is None:
            self.isoc_args = {}
        if make_isocs:
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


class ExtendedSolarIsocs(IsochroneSetBase):
    """Solar metallicity isochrone set."""
    def __init__(self, **kwargs):
        self.isoc_args = kwargs.pop('isoc_args', dict())
        self.isoc_phases = kwargs.pop('isoc_phases', None)
        self.z_grid = (0.0096, 0.012, 0.015, 0.019, 0.024, 0.030, 0.04)
        print "ExtendedSolarZIsocs", kwargs
        super(ExtendedSolarIsocs, self).__init__(**kwargs)

    @property
    def bands(self):
        return ('F475W', 'F814W', 'F275W', 'F336W', 'F110W', 'F160W')

    @property
    def distance(self):
        return Distance(785. * u.kpc)

    def setup_isochrones(self):
        """Download Padova isochrones."""
        print "Running ExtendedSolarIsocs setup_isochrones"
        WFC3_BANDS = ['F275W1', 'F336W', 'F110W', 'F160W']
        ACS_BANDS = ['F475W', 'F814W']

        if not os.path.exists(os.path.join(STARFISH, self.isoc_dir)):
            make_isocs = True
        elif len(glob(os.path.join(STARFISH, self.isoc_dir, "*"))) == 0:
            make_isocs = True
        else:
            make_isocs = False

        if self.isoc_args is None:
            self.isoc_args = {}
        if make_isocs:
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
        super(ExtendedSolarIsocs, self).setup_isochrones()


class PhatCatalog(DatasetBase):
    """Mixin for PHAT photometry data.

    Photometry is lazy loaded to it is efficient to rebuild the pipeline
    object.
    """
    def __init__(self, brick, **kwargs):
        self.brick = brick
        self._phat_data = None
        super(PhatCatalog, self).__init__(**kwargs)

    def _load_phat_data(self):
        self._phat_data = Table.read(phat_v2_phot_path(self.brick),
                                     format='fits')
        phat_bands = ('F475W', 'F814W', 'F275W', 'F336W', 'F110W', 'F160W')
        # Normalize bandpass names
        for band in phat_bands:
            old_name = "_".join((band.lower(), 'vega'))
            self._phat_data.rename_column(old_name, band)

    def get_phot(self, band):
        if self._phat_data is None:
            self._load_phat_data()  # lazy loading

        if not isinstance(band, basestring):
            band1, band2 = band
            return self._phat_data[band1] - self._phat_data[band2]
        else:
            return self._phat_data[band]

    def write_phot(self, x_band, y_band, data_root, suffix):
        if self._phat_data is None:
            self._load_phat_data()  # lazy loading

        x = self.get_phot(x_band)
        y = self.get_phot(y_band)
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


class ExtendedSolarLockfile(LockBase):
    """Lockfile mixin to create an sub-, solar and super solar Z Hess set."""
    def __init__(self, **kwargs):
        print "ExtendedSolarLockfile", kwargs
        super(ExtendedSolarLockfile, self).__init__(**kwargs)

    def build_lockfile(self):
        self.lockfile = Lockfile(self.builder.read_isofile(), self.synth_dir,
                                 unbinned=False)

        z_bins = [(0.009, 0.0135), (0.0135, 0.022), (0.022, 0.032)]
        z_strs = ['0010', '0019', '0028']

        # Bin young isochrones
        for zbin, z_str in zip(z_bins, z_strs):
            young_grid = np.linspace(6.5, 8.95, 10)
            for i, logage0 in enumerate(young_grid[:-1]):
                logage0 = logage0
                logage1 = young_grid[i + 1]
                mean_age = (logage0 + logage1) / 0.2
                name = "z{0}_{1:05.2f}".format(z_str, mean_age)
                self.lockfile.lock_box(name, (logage0, logage1), zbin)

            # Bin old isochrones
            old_grid = np.arange(1e9, 14 * 1e9, 1e9)
            for i, age0 in enumerate(old_grid[:-1]):
                logage0 = np.log10(age0 - 0.05 * 1e9)
                logage1 = np.log10(old_grid[i + 1])
                mean_age = (logage0 + logage1) / 0.2
                name = "z{0}_{1:05.2f}".format(z_str, mean_age)
                self.lockfile.lock_box(name, (logage0, logage1), zbin)


class PhatCrowding(CrowdingBase):
    """Use crowding from the PHAT AST fields."""
    def __init__(self, **kwargs):
        self._ast_field = kwargs.pop('ast_field', 0)
        super(PhatCrowding, self).__init__(**kwargs)

    def build_crowding(self):
        # Use PHAT AST from the outer field (field 0)
        crowd_path = os.path.join(self.synth_dir, "crowding.dat")
        full_crowd_path = os.path.join(STARFISH, crowd_path)
        tbl = PhatAstTable()
        tbl.write_crowdfile_for_field(full_crowd_path, 0,
                                      bands=self.bands)
        self.crowd = ExtantCrowdingTable(crowd_path)


class NoDust(ExtinctionBase):
    """Mixin for no dust."""
    def __init__(self, **kwargs):
        super(NoDust, self).__init__(**kwargs)

    def build_extinction(self):
        # Add MW dust screen
        self.young_av = ExtinctionDistribution()
        self.young_av.set_uniform(mw_Av())
        self.old_av = ExtinctionDistribution()
        self.old_av.set_uniform(mw_Av())

        self.rel_extinction = phat_rel_extinction()


class PhatGaussianDust(ExtinctionBase):
    """Mixin for Gaussian dust distributions for PHAT filters."""
    def __init__(self, **kwargs):
        self._young_av = kwargs.pop('young_av', 1.0)
        self._old_av = kwargs.pop('old_av', 0.5)
        self._av_sigma_ratio = kwargs.pop('av_sigma_ratio', 0.5)
        super(PhatGaussianDust, self).__init__(**kwargs)

    def build_extinction(self):
        # NOTE includes the E(V
        self.young_av = ExtinctionDistribution()
        if self._young_av > 0.:
            av = np.random.normal(
                loc=self._young_av,
                scale=self._young_av * self._av_sigma_ratio,
                size=1000)
            av[av < 0.] = 0.
            self.young_av.set_samples(av + mw_Av())
        else:
            self.young_av.set_samples(np.zeros(1000) + mw_Av())

        self.old_av = ExtinctionDistribution()
        if self._old_av > 0.:
            av = np.random.normal(
                loc=self._old_av,
                scale=self._old_av * self._av_sigma_ratio,
                size=1000)
            av[av < 0.] = 0.
            self.old_av.set_samples(av + mw_Av())
        else:
            self.old_av.set_samples(np.zeros(1000) + mw_Av())

        self.rel_extinction = phat_rel_extinction()


class PhatStepDust(ExtinctionBase):
    """Mixin for Uniform dust distributions for PHAT filters."""
    def __init__(self, **kwargs):
        self._young_av = kwargs.pop('young_av', 0.0) + mw_Av()
        self._old_av = kwargs.pop('old_av', 0.0) + mw_Av()
        self._young_dav = kwargs.pop('young_dav', 1.0)
        self._old_dav = kwargs.pop('old_dav', 1.0)
        super(PhatStepDust, self).__init__(**kwargs)

    def build_extinction(self):
        # NOTE includes MW extinction from __init__
        self.young_av = ExtinctionDistribution()
        self.young_av.set_samples(np.random.uniform(
            low=self._young_av,
            high=self._young_av + self._young_dav,
            size=1000))

        self.old_av = ExtinctionDistribution()
        self.old_av.set_samples(np.random.uniform(
            low=self._old_av,
            high=self._old_av + self._old_dav,
            size=1000))

        self.rel_extinction = phat_rel_extinction()


class LewisBrickDust(ExtinctionBase):
    """Mixin for a uniform dust distribution fitted from Lewis et al 15 Fig 17.

    The maximum extinction is estimated from a Draine et al 2015 dust map.
    Requires that the brick be known.
    """
    def __init__(self):
        super(LewisBrickDust, self).__init__()

    def build_extinction(self):
        """Young and old dust at equal here."""
        # Get the coordinate of the brick
        brick_fits = phat_brick_path(self.brick, 'F814W')
        wcs = astropy.io.WCS(brick_fits[0].header)
        poly = wcs.calc_footprint()
        ra0 = poly[0, :].mean()
        dec0 = poly[1, :].mean()

        lewis = LewisDustLaw()
        max_av = lewis.estimate_extinction(ra0, dec0)
        av = np.random.uniform(low=mw_Av(),
                               high=max_av,
                               size=1000)

        self.young_av = ExtinctionDistribution()
        self.young_av.set_samples(av)

        self.old_av = ExtinctionDistribution()
        self.old_av.set_samples(av)

        self.rel_extinction = phat_rel_extinction()


class SolarZPhatPipeline(CompletePhatPlanes, SolarZIsocs,
                         SolarLockfile, NoDust, PhatCrowding, PipelineBase):
    """A pipeline for fitting PHAT bricks with solar metallicity isochrones."""
    def __init__(self, **kwargs):
        print "SolarZPhatPipeline", kwargs
        super(SolarZPhatPipeline, self).__init__(**kwargs)


class SolarRgbPipeline(RgbPhatPlanes, SolarZIsocs,
                       SolarLockfile, NoDust, PhatCrowding, PipelineBase):
    """A pipeline for fitting PHAT bricks with solar metallicity isochrones."""
    def __init__(self, **kwargs):
        print "SolarRgbPipeline", kwargs
        super(SolarRgbPipeline, self).__init__(**kwargs)


def build_phat_filter_set(**kwargs):
    r_wfc3 = AgeGridRequest(photsys='wfc3_wide', **kwargs)
    r_acs = AgeGridRequest(photsys='acs_wfc', **kwargs)
    isoc_set = join_isochrone_sets(r_wfc3.isochrone_set,
                                   r_acs.isochrone_set,
                                   left_bands=['F275W1', 'F336W',
                                               'F110W', 'F160W'],
                                   right_bands=['F475W', 'F814W'])
    return isoc_set


get_demo_age_grid = partial(build_phat_filter_set,
                            z=0.019, min_log_age=6.6, max_log_age=10.13,
                            delta_log_age=0.2)


def plot_isochrone_phases(ax, band1, band2, show_cb=False, cb_ax=None):
    isoc_set = get_demo_age_grid(**dict(isoc_kind='parsec_CAF09_v1.2S',
                                        photsys_version='yang'))
    phase_labels = {0: 'Pre-MS', 1: 'MS', 2: 'SGB', 3: 'RGB',
                    4: 'CHeB(1)', 5: 'CHeB(2)', 6: 'CHeB(3)',
                    7: 'E-AGB', 8: 'TP-AGB'}
    cmap = mpl.colors.ListedColormap(
        palettable.colorbrewer.qualitative.Set1_9.mpl_colors)
    scalar_map = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=-0.5,
                                                                 vmax=8.5),
                                       cmap=cmap)
    scalar_map.set_array(np.array(range(0, 9)))

    d = Distance(785 * u.kpc)
    for isoc in isoc_set:
        phases = np.unique(isoc['stage'])
        srt = np.argsort(phases)
        phases = phases[srt]
        for p in phases:
            s = np.where(isoc['stage'] == p)[0]
            ax.plot(isoc[band1][s] - isoc[band2][s],
                    isoc[band2][s] + d.distmod.value,
                    c=scalar_map.to_rgba(p),
                    lw=0.8)
    if show_cb:
        cb = plt.colorbar(mappable=scalar_map,
                          cax=cb_ax, ax=ax, ticks=range(0, 9))
        cb.ax.set_yticklabels([phase_labels[p] for p in range(0, 9)])
        cb.set_label(r"Stage")
