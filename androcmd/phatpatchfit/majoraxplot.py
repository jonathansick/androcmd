# encoding: utf-8
"""
Tools for plotting SFR along the major axis.

2015-07-01 - Created by Jonathan Sick
"""

import matplotlib as mpl
import numpy as np

from .analysistools import marginalize_metallicity


def select_patches(dataset):
    """Patch keys that lie along the major axis."""
    keys = []
    phi = []
    for key, patch in dataset['patches'].items():
        keys.append(key)
        phi.append(patch.attrs['phi'])
    phi = np.array(phi)
    sel = np.where((phi <= 20.) | (phi >= 340.))[0]
    return [keys[i] for i in sel]


def bin_patches_radially(dataset, patch_keys):
    r_grid = np.arange(0, 21, 1.)
    binned_patches = [[] for r in r_grid]
    patch_r = np.array([dataset['patches'][k].attrs['r_kpc']
                        for k in patch_keys])
    for i in xrange(len(r_grid) - 1):
        rmin = r_grid[i]
        rmax = r_grid[i + 1]
        s = np.where((patch_r >= rmin) & (patch_r < rmax))[0]
        for si in s:
            binned_patches[i].append(patch_keys[si])
    return r_grid, binned_patches


def plot_highlighted_patches(dataset, patch_keys, ax):
    for key in patch_keys:
        patch = dataset['patches'][key]
        poly = patch.attrs['poly']
        patch = mpl.patches.Polygon(poly, closed=True,
                                    transform=ax.get_transform('world'),
                                    facecolor='y', alpha=0.5,
                                    edgecolor='k', lw=0.5)
        ax.add_patch(patch)


def compute_sfr_in_span(dataset, patch_keys, fit_key, myr_min, myr_max):
    """Compute the mean SFR of all patches in a bin within a time span (in Myr)
    """
    patch_sfrs = []
    for k in patch_keys:
        logage, sfr = marginalize_metallicity(dataset['patches'][k], fit_key)
        a_myr = 10. ** (logage - 6)
        interp_ages = np.linspace(myr_min, myr_max, 50)
        interp_sfrs = np.interp(interp_ages, a_myr, sfr)
        mean_sfr = np.mean(interp_sfrs)
        patch_sfrs.append(mean_sfr)
    patch_sfrs = np.array(patch_sfrs)
    if len(patch_sfrs) == 0:
        return None
    else:
        mean_sfr = patch_sfrs.mean()
        return mean_sfr
