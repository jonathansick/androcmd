#!/usr/bin/env python
# encoding: utf-8
"""
Tools for working with dust attenuation for M31.

2015-05-12 - Created by Jonathan Sick
"""

import numpy as np
import astropy

from m31hst.draine import spire350_dust_mass_map
from starfisher.dust import SF11ExtinctionCurve


def mw_Av():
    """Build the A_V attenuation by the MW towards M31."""
    curve = SF11ExtinctionCurve()
    ratio = curve['Landolt V']  # A_V / E(B-V) from T6 of SF2011
    return ratio * 0.07  # 0.07 is E(B-V) to M31 from Schlegel 1998


def phat_rel_extinction():
    """Build array of relative extinction in PHAT bands."""
    curve = SF11ExtinctionCurve()
    rel_extinction = np.ones(6, dtype=float)
    rel_extinction[0] = curve.extinction_ratio('WFC3 F275W')
    rel_extinction[1] = curve.extinction_ratio('WFC3 F336W')
    rel_extinction[2] = curve.extinction_ratio('ACS F475W')
    rel_extinction[3] = curve.extinction_ratio('ACS F814W')
    rel_extinction[4] = curve.extinction_ratio('WFC3 F110W')
    rel_extinction[5] = curve.extinction_ratio('WFC3 F160W')
    return rel_extinction


class LewisDustLaw(object):
    """Assign extinction Av+dAv from the Drain et al dust map."""
    def __init__(self):
        super(LewisDustLaw, self).__init__()
        self._f = astropy.io.fits.open(spire350_dust_mass_map())
        self._wcs = astropy.wcs.WCS(self._f[0].header)
        self._naxis1 = self._f[0].header['NAXIS1']
        self._naxis2 = self._f[0].header['NAXIS2']

    def estimate_extinction(self, ra, dec):
        """Estimate the maximum Av estinction to this coordinate."""
        x, y = self.wcs.all_world2pix([ra], [dec], 0)
        assert x > 0 and x < self._naxis1
        assert y > 0 and y < self._naxis2
        x = int(x)
        y = int(y)
        Sigma_dust = self._f[0].data[y, x]
        return 10. ** (-5.4) * Sigma_dust
