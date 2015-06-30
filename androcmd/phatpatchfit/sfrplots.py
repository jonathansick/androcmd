# encoding: utf-8
"""
Tools for plotting SFR

2015-06-30 - Created by Jonathan Sick
"""

import numpy as np

from .analysistools import marginalize_metallicity


SFR_LABEL = r'$\langle \mathrm{SFR} \rangle$ ($\mathrm{M}_\odot' \
            r'~10^{-3}~\mathrm{yr}^{-1}~\mathrm{kpc}^{-2}$)'


def get_scaled_sfr_values(dataset, fit_key, age):
    patches = dataset['patches']

    ra = []
    dec = []
    log_sfrs = []
    for patch_name, patch_group in patches.items():
        logage_tbl, sfr_tbl = marginalize_metallicity(patch_group, fit_key)
        sfr = np.interp(np.log10(age), logage_tbl, sfr_tbl)
        log_sfrs.append(scale_sfr(sfr, patch_group))
        ra.append(patch_group.attrs['ra0'])
        dec.append(patch_group.attrs['dec0'])

    return np.array(ra), np.array(dec), np.array(log_sfrs)


def scale_sfr(sfr, patch_group):
    area = patch_group.attrs['area_proj'] \
        / np.cos(77.5 * np.pi / 180.) / 1e3 / 1e3  # kpc^2
    return np.log10(sfr / area * 10. ** 3.)
