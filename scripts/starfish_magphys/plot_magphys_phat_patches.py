#!/usr/bin/env python
# encoding: utf-8
"""
Plot MAGPHYS estimates of PHAT fields.

2015-07-16 - Created by Jonathan Sick
"""

import argparse
import os

import h5py
import matplotlib as mpl
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.gridspec as gridspec
from palettable.cubehelix import perceptual_rainbow_16

from androcmd.phatpatchfit.pipeline import load_field_patches
from androcmd.phatpatchfit import load_galex_map, setup_galex_axes


def main():
    args = parse_args()
    base_path = os.path.splitext(args.path)[0]

    plot_mean_age(args.path, base_path + '_mean_age')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('path',
                        help='andromass pixel table with MAGPHYS estimates')
    return parser.parse_args()


def plot_mean_age(sed_path, plot_path):
    basemap = load_galex_map()

    dataset = h5py.File(sed_path, 'r')
    t = dataset['estimates']
    log_age_M = t['log_age_M']
    age_gyr = 10. ** (log_age_M - 9.)
    polys = _make_footprint_list()

    fig = Figure(figsize=(3.5, 3.5), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(1, 2,
                           left=0.12, right=0.85, bottom=0.15, top=0.95,
                           wspace=0., hspace=None,
                           width_ratios=(1, 0.08), height_ratios=None)
    ax = setup_galex_axes(fig, gs[0], basemap)
    ax_cb = fig.add_subplot(gs[1])
    print age_gyr[:, 2].flatten()
    cmap = perceptual_rainbow_16.mpl_colormap
    norm = mpl.colors.Normalize(vmin=0, vmax=10, clip=True)
    mapper = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    mapper.set_array(age_gyr[:, 2].flatten())

    poly_args = {'edgecolor': 'k', 'lw': 0.5,
                 'transform': ax.get_transform('world'),
                 'closed': True}
    for poly, age in zip(polys, age_gyr):
        patch = mpl.patches.Polygon(poly,
                                    facecolor=mapper.to_rgba(age[2]),
                                    **poly_args)
        ax.add_patch(patch)
    cb = fig.colorbar(mapper, cax=ax_cb)
    cb.set_label(r'$\langle A~\mathrm{Gyr}^{-1} \rangle$')

    ax.coords[0].ticklabels.set_size(8)
    ax.coords[1].ticklabels.set_size(8)

    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def _make_footprint_list():
    """Make a list of footprints that matches the SED ID order"""
    fields = load_field_patches()
    return [patch['poly'] for patch in fields]


if __name__ == '__main__':
    main()
