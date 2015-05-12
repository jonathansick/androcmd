#!/usr/bin/env python
# encoding: utf-8
"""
Tools for working with dust attenuation for M31.

2015-05-12 - Created by Jonathan Sick
"""

import numpy as np

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
