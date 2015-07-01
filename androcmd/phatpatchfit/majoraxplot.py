# encoding: utf-8
"""
Tools for plotting SFR along the major axis.

2015-07-01 - Created by Jonathan Sick
"""

import matplotlib as mpl
import numpy as np


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


def plot_highlighted_patches(dataset, patch_keys, ax):
    for key in patch_keys:
        patch = dataset['patches'][key]
        poly = patch.attrs['poly']
        patch = mpl.patches.Polygon(poly, closed=True,
                                    transform=ax.get_transform('world'),
                                    facecolor='y', alpha=0.5,
                                    edgecolor='k', lw=0.5)
        ax.add_patch(patch)
