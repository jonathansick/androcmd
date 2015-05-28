# encoding: utf-8
"""
Tools for the pan-M31 PHAT patch fitting exercise.

2015-05-28 - Created by Jonathan Sick
"""

import os
from collections import OrderedDict

from matplotlib.patches import Polygon

import numpy as np
from astropy.table import Table

from m31hst import phat_v2_phot_path

from starfisher.pipeline import PipelineBase
from starfisher.pipeline import PlaneBase
from starfisher.pipeline import ExtinctionBase
from starfisher.pipeline import DatasetBase
from starfisher.dust import ExtinctionDistribution

# from starfisher.pipeline import (
#     PipelineBase, IsochroneSetBase, DatasetBase, LockBase,
#     CrowdingBase, ExtinctionBase)

from androcmd.planes import make_f475w_f160w

from androcmd.phatpipeline import (PhatCrowding,
                                   ExtendedSolarIsocs, ExtendedSolarLockfile)
from androcmd.dust import mw_Av, phat_rel_extinction, LewisDustLaw


class BaselineTestPhatPlanes(PlaneBase):
    """Color plane set for the PHAT color baseline comparison test."""
    def __init__(self, **kwargs):
        self._planes = OrderedDict([
            ('oir_all', make_f475w_f160w()),
        ])
        super(BaselineTestPhatPlanes, self).__init__(**kwargs)

    @property
    def planes(self):
        return self._planes


class LewisPatchDust(ExtinctionBase):
    """Mixin for a uniform dust distribution fitted from Lewis et al 15 Fig 17.

    The maximum extinction is estimated from a Draine et al 2015 dust map.
    Requires that the brick be known.
    """
    def __init__(self, **kwargs):
        # self.brick = kwargs.pop('brick', 23)
        super(LewisPatchDust, self).__init__(**kwargs)

    def build_extinction(self):
        """Young and old dust at equal here."""
        lewis = LewisDustLaw()
        max_av = lewis.estimate_mean_extinction(self.poly)
        av = np.random.uniform(low=mw_Av(),
                               high=max_av,
                               size=1000)

        self.young_av = ExtinctionDistribution()
        self.young_av.set_samples(av)

        self.old_av = ExtinctionDistribution()
        self.old_av.set_samples(av)

        self.rel_extinction = phat_rel_extinction()


class ThreeZPipeline(BaselineTestPhatPlanes, ExtendedSolarIsocs,
                     ExtendedSolarLockfile, LewisPatchDust, PhatCrowding,
                     PipelineBase):
    """Pipeline for patch fitting with three metallicity tracks."""
    def __init__(self, **kwargs):
        # Get patch attributes
        self.patch = kwargs.pop('patch')
        self.poly = kwargs.pop('poly')
        self.brick = kwargs.pop('brick')
        self.ra0 = kwargs.pop('ra0')
        self.dec0 = kwargs.pop('dec0')
        self.area = kwargs.pop('area')
        super(ThreeZPipeline, self).__init__(**kwargs)


class PatchCatalog(DatasetBase):
    """Mixin for PHAT photometry data.

    Photometry is lazy loaded to it is efficient to rebuild the pipeline
    object.
    """
    def __init__(self, **kwargs):
        self.patch = kwargs.pop('patch')
        self.poly = kwargs.pop('poly')
        self.brick = kwargs.pop('brick')
        self.ra0 = kwargs.pop('ra0')
        self.dec0 = kwargs.pop('dec0')
        self.area = kwargs.pop('area')
        super(PatchCatalog, self).__init__(**kwargs)

    def _load_phat_data(self):
        self._phat_data = Table.read(phat_v2_phot_path(self.brick),
                                     format='fits')
        phat_bands = ('F475W', 'F814W', 'F275W', 'F336W', 'F110W', 'F160W')
        # Normalize bandpass names
        for band in phat_bands:
            old_name = "_".join((band.lower(), 'vega'))
            self._phat_data.rename_column(old_name, band)

        # Spatially select stars in the patch footprint
        self._phat_data = self._spatial_selection(self._phat_data)

    def _spatial_selection(self, dataset):
        data_ra = dataset['ra']
        data_dec = dataset['dec']
        xy = np.vstack((data_ra, data_dec)).T
        print "xy.shape", xy.shape

        polygon = Polygon(self.poly, closed=True)
        inside = polygon.contains_point(xy)
        return dataset[inside == True]  # noqa

    def get_phot(self, band):
        if self._phat_data is None:
            self._load_phat_data()  # lazy loading

        if not isinstance(band, basestring):
            band1, band2 = band
            return self._phat_data[band1] - self._phat_data[band2]
        else:
            return self._phat_data[band]

    def _select_gst(self, x_band, y_band):
        mags = []
        if not isinstance(x_band, basestring):
            mags.extend(x_band)
        else:
            mags.append(x_band)

        if not isinstance(y_band, basestring):
            mags.extend(y_band)
        else:
            mags.append(y_band)
        mags = list(set(mags))

        gsts = []
        for band in mags:
            key = '{0}_gst'.format(band.lower())
            gsts.append(self._phat_data[key])
        gsts_array = np.vstack(gsts).T
        gsts = np.all(gsts_array, axis=1)
        return np.where(gsts == True)[0]  # NOQA

    def write_phot(self, x_band, y_band, data_root, suffix):
        """Only good (GST=1 in all relevant bands) photometry is written."""
        if self._phat_data is None:
            self._load_phat_data()  # lazy loading

        x = self.get_phot(x_band)
        y = self.get_phot(y_band)
        gst_sel = self._select_gst(x_band, y_band)

        phot_dtype = np.dtype([('x', np.float), ('y', np.float)])
        photdata = np.empty(len(gst_sel), dtype=phot_dtype)
        photdata['x'][:] = x[gst_sel]
        photdata['y'][:] = y[gst_sel]

        path = data_root + suffix
        full_path = os.path.join(os.getenv('STARFISH'), path)
        fit_dir = os.path.dirname(full_path)
        if not os.path.exists(fit_dir):
            os.makedirs(fit_dir)
        np.savetxt(full_path, photdata, delimiter=' ', fmt='%.4f')

    @property
    def polygon(self):
        """Polygon bounding box around the dataset."""
        return self.poly
