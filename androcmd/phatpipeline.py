#!/usr/bin/env python
# encoding: utf-8
"""
Pipelines for fitting PHAT CMDs.
"""
import os
from glob import glob

import numpy as np
from astropy.coordinates import Distance
import astropy.units as u
from astropy.table import Table

from m31hst import phat_v2_phot_path
from m31hst.phatast import PhatAstTable

from padova import AgeGridRequest
from padova.isocdata import join_isochrone_sets, Isochrone

from starfisher import Lockfile
from starfisher import ExtantCrowdingTable
from starfisher.pipeline import (
    PipelineBase, IsochroneSetBase, DatasetBase, LockBase,
    CrowdingBase, ExtinctionBase)

from androcmd.planes import RgbPhatPlanes, CompletePhatPlanes


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
            old_name = "_".join((band.lower(), 'vega'))
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
        super(NoDust, self).build_extinction()


class SolarZPhatPipeline(CompletePhatPlanes, SolarZIsocs, PhatCatalog,
                         SolarLockfile, NoDust, PhatCrowding, PipelineBase):
    """A pipeline for fitting PHAT bricks with solar metallicity isochrones."""
    def __init__(self, **kwargs):
        print "SolarZPhatPipeline", kwargs
        super(SolarZPhatPipeline, self).__init__(**kwargs)


class SolarRgbPipeline(RgbPhatPlanes, SolarZIsocs, PhatCatalog,
                       SolarLockfile, NoDust, PhatCrowding, PipelineBase):
    """A pipeline for fitting PHAT bricks with solar metallicity isochrones."""
    def __init__(self, **kwargs):
        print "SolarRgbPipeline", kwargs
        super(SolarRgbPipeline, self).__init__(**kwargs)
