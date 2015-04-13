#!/usr/bin/env python
# encoding: utf-8
"""
Pipeline for analysis of hte PHAT six-color photometry.
"""

import os
import numpy as np
from astropy.table import Table
from m31hst import phat_v2_phot_path


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
