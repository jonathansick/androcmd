# encoding: utf-8
"""
Tools for plotting SFR

2015-06-30 - Created by Jonathan Sick
"""

import numpy as np


SFR_LABEL = r'$\langle \log_{10} \Sigma_\mathrm{SFR} \rangle ~'\
            r'(10^{-3}~\mathrm{M}_\odot~\mathrm{yr}^{-1}'\
            r'~\mathrm{kpc}^{-2}$)'


LIN_SFR_LABEL = r'$\langle \Sigma_\mathrm{SFR} \rangle ~'\
                r'(10^{-3}~\mathrm{M}_\odot~\mathrm{yr}^{-1}'\
                r'~\mathrm{kpc}^{-2}$)'


def get_scaled_sfr_values(dataset, fit_key, age):
    patches = dataset['patches']

    ra = []
    dec = []
    log_sfrs = []
    for patch_name, patch_group in patches.items():
        t = patch_group['sfh_marginal'][fit_key]
        logage_tbl = t['log(age)']
        sfr_tbl = t['sfr']
        sfr = np.interp(np.log10(age), logage_tbl, sfr_tbl)
        log_sfrs.append(scale_sfr(sfr, patch_group))
        ra.append(patch_group.attrs['ra0'])
        dec.append(patch_group.attrs['dec0'])

    return np.array(ra), np.array(dec), np.array(log_sfrs)


def scale_sfr(sfr, patch_group):
    area = patch_group.attrs['area_proj'] \
        / np.cos(77.5 * np.pi / 180.) / 1e3 / 1e3  # kpc^2
    return np.log10(sfr / area * 10. ** 3.)


def lin_scale_sfr(sfr, patch_group):
    area = patch_group.attrs['area_proj'] \
        / np.cos(77.5 * np.pi / 180.) / 1e3 / 1e3  # kpc^2
    return sfr / area * 10. ** 3.
