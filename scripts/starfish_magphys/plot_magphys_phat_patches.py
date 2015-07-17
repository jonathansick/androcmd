#!/usr/bin/env python
# encoding: utf-8
"""
Plot MAGPHYS estimates of PHAT fields.

2015-07-16 - Created by Jonathan Sick
"""

import argparse
import os

import numpy as np
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

    # plot_mean_age(args.path, base_path + '_mean_age')
    plot_estimates_grid(args.path, base_path + '_est_grid')


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


def plot_estimates_grid(sed_path, plot_path):
    dataset = h5py.File(sed_path, 'r')

    fig = Figure(figsize=(6.5, 7), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(2, 3,
                           left=0.04, right=0.96, bottom=0.07, top=0.95,
                           wspace=0.15, hspace=0.15,
                           width_ratios=None, height_ratios=None)
    gs_settings = dict(height_ratios=(1, 0.05), hspace=0.01)
    gs_age = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[0, 0],
                                              **gs_settings)
    gs_sfr = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1, 0],
                                              **gs_settings)
    gs_Z = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[0, 1],
                                            **gs_settings)
    gs_mu = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1, 1],
                                             **gs_settings)
    gs_tau = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[0, 2],
                                              **gs_settings)
    gs_ism = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1, 2],
                                              **gs_settings)

    t = dataset['estimates']
    _plot_grid_ax(10. ** (t['log_age_M'] - 9.), fig, gs_age,
                  vmin=0., vmax=12.,
                  label=r'$\langle A~\mathrm{Gyr}^{-1} \rangle$')

    _plot_grid_ax(t['Z_Zo'], fig, gs_Z,
                  vmin=-2, vmax=0.3,
                  label=r'$\log Z/Z_\odot$')

    print t.dtype
    _plot_grid_ax(scale_sfr(t['SFR_0.1Gyr']), fig, gs_sfr,
                  vmin=1.0, vmax=1.25,
                  label=r'$\langle \log \mathrm{SFR}_{100~\mathrm{Myr}} \rangle$')

    _plot_grid_ax(t['mu'], fig, gs_mu,
                  vmin=0., vmax=1,
                  label=r'$\mu$')

    _plot_grid_ax(t['tau_V'], fig, gs_tau,
                  vmin=0., vmax=2,
                  label=r'$\tau_V$')

    _plot_grid_ax(t['tau_V_ISM'], fig, gs_ism,
                  vmin=0., vmax=0.5,
                  label=r'$\tau_\mathrm{ISM}$')

    canvas.print_figure(plot_path + ".pdf", format="pdf")


def _plot_grid_ax(data, fig, gs, vmin=0., vmax=1.,
                  label=None, cb_locator=mpl.ticker.MaxNLocator(4)):
    basemap = load_galex_map()
    polys = _make_footprint_list()

    ax = setup_galex_axes(fig, gs[0], basemap)
    ax_cb = fig.add_subplot(gs[1])
    print data[:, 2].flatten()
    cmap = perceptual_rainbow_16.mpl_colormap
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
    mapper = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    mapper.set_array(data[:, 2].flatten())

    poly_args = {'edgecolor': 'k', 'lw': 0.5,
                 'transform': ax.get_transform('world'),
                 'closed': True}
    for poly, v in zip(polys, data):
        patch = mpl.patches.Polygon(poly,
                                    facecolor=mapper.to_rgba(v[2]),
                                    **poly_args)
        ax.add_patch(patch)
    cb = fig.colorbar(mapper, cax=ax_cb, orientation='horizontal')
    if label is not None:
        cb.set_label(label, fontsize=9)
        cb.ax.tick_params(labelsize=7)
        cb.locator = cb_locator
        cb.update_ticks()

    ax.coords[0].ticklabels.set_visible(False)
    ax.coords[1].ticklabels.set_visible(False)


def _make_footprint_list():
    """Make a list of footprints that matches the SED ID order"""
    fields = load_field_patches()
    return [patch['poly'] for patch in fields]


def scale_sfr(sfr):
    fields = load_field_patches()
    area_proj = np.atleast_2d(np.array([f['area_proj'] for f in fields])).T
    print "area_proj"
    print area_proj
    area = area_proj / np.cos(77.5 * np.pi / 180.) / 1e3 / 1e3  # kpc^2
    print "area"
    print area
    print "sfr"
    print sfr
    scaled_sfr = np.log10(10. ** sfr / area * 10. ** 3.)
    print "scaled log(sfr)"
    print scaled_sfr
    return scaled_sfr


if __name__ == '__main__':
    main()
