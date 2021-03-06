#!/usr/bin/env python
# encoding: utf-8
"""
Plotting and analysis of the PHAT patch fitting.

Use the fetch_patch_results.py script to build a patch dataset.

2015-06-17 - Created by Jonathan Sick
"""

import argparse
import h5py
import numpy as np

import matplotlib as mpl
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.gridspec as gridspec
import palettable
from palettable.cubehelix import perceptual_rainbow_16

from androcmd.phatpatchfit import (load_galex_map,
                                   setup_galex_axes, setup_plane_comp_axes,
                                   setup_plane_axes,
                                   get_scaled_sfr_values,
                                   scale_sfr,
                                   SFR_LABEL, LIN_SFR_LABEL,
                                   plot_patch_footprints)
import androcmd.phatpatchfit.majoraxplot as majoraxplot


def main():
    args = parse_args()

    dataset = h5py.File(args.hdf5_path, 'r')

    if args.radial_mean_age is not None:
        plot_radial_mean_age(dataset, args.radial_mean_age)

    if args.rchi_hist is not None:
        plot_reduced_chi_hist(dataset, args.rchi_hist)

    if args.sfh_lines is not None:
        plot_sfh_lines(dataset, args.sfh_lines)

    if args.mean_age_map is not None:
        plot_mean_age_map(dataset, args.mean_age_map)
        plot_single_mean_age_map(dataset, 'oir_all',
                                 args.mean_age_map + '_oir_all')

    if args.epoch_sfr_maps is not None:
        for fit_key in ['lewis', 'oir_all']:
            plot_epoch_sfr_map_vertical(
                dataset, fit_key,
                '_'.join((args.epoch_sfr_maps, fit_key, 'vert')))

    if args.major_ax_sfr is not None:
        plot_major_ax_sfr(dataset, args.major_ax_sfr)
        plot_major_ax_sfr_linear(dataset, args.major_ax_sfr + '_linear')

    if args.major_ax_cmass is not None:
        plot_major_ax_cumulative_mass_unnorm(
            dataset, args.major_ax_cmass + '_unnorm')
        plot_major_ax_cumulative_mass_norm(
            dataset, args.major_ax_cmass + '_norm')

    if args.cmass_age_map is not None:
        plot_mass_accumulation_age_map(
            dataset,
            args.cmass_age_map,
            fit_key='oir_all')


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('hdf5_path')
    parser.add_argument('--radial-mean-age', default=None)
    parser.add_argument('--rchi-hist', default=None)
    parser.add_argument('--sfh-lines', default=None)
    parser.add_argument('--mean-age-map', default=None)
    parser.add_argument(
        '--epoch-sfr-maps', default=None,
        help='Map comparing SFR at for MS and ALL, individually')
    parser.add_argument('--major-ax-sfr', default=None,
                        help='Emulate Lewis 2015 Fig 6')
    parser.add_argument(
        '--major-ax-cmass', default=None,
        help='Major axis cumulative mass build-up')
    parser.add_argument(
        '--cmass-age-map', default=None,
        help='Map of ages when cumulative mass thresholds are reached')
    return parser.parse_args()


def plot_radial_mean_age(dataset, plot_path):
    fig = Figure(figsize=(3.5, 3.5), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(1, 1,
                           left=0.15, right=0.95, bottom=0.15, top=0.95,
                           wspace=None, hspace=None,
                           width_ratios=None, height_ratios=None)
    ax = fig.add_subplot(gs[0])
    R = dataset['sfh_table']['r_kpc'][:]
    mean_age = dataset['sfh_table']['mean_age_oir_all'][:]
    # mean_age_ms = dataset['sfh_table']['mean_age_lewis'][:]
    ax.scatter(R, mean_age, s=3, edgecolors='None', facecolors='k',
               label='OIR-ALL')
    # ax.scatter(R, mean_age_ms, s=3, edgecolors='None', facecolors='dodgerblue',
    #            label='ACS-MS')
    ax.set_xlim(0., 25.)
    ax.set_ylim(0., 10.)
    ax.set_xlabel(r'$R_\mathrm{maj}$ (kpc)')
    ax.set_ylabel(r'$\langle A \rangle$ (Gyr)')
    # ax.legend(frameon=True, markerscale=2, ncol=2)
    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def plot_reduced_chi_hist(dataset, plot_path):
    chisq_ms = dataset['sfh_table']['chi_red_lewis'][:]
    chisq_oir = dataset['sfh_table']['chi_red_oir_all'][:]

    fig = Figure(figsize=(3.5, 3.5), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(1, 1,
                           left=0.15, right=0.95, bottom=0.15, top=0.95,
                           wspace=None, hspace=None,
                           width_ratios=None, height_ratios=None)
    ax = fig.add_subplot(gs[0])
    ax.hist(chisq_ms, 20, histtype='step', edgecolor='dodgerblue',
            label='ACS-MS')
    ax.hist(chisq_oir, 20, histtype='step', edgecolor='firebrick',
            label='OIR-ALL')
    ax.set_xlabel(r'$\chi_r^2$')
    ax.set_ylabel(r'Patches')
    ax.set_xlim(0, 30)
    ax.legend()
    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def plot_sfh_lines(dataset, plot_path):
    fig = Figure(figsize=(6, 6), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(2, 2,
                           left=0.1, right=0.89, bottom=0.1, top=0.95,
                           wspace=0.1, hspace=0.1,
                           width_ratios=(1, 0.1), height_ratios=None)
    ax_ms = fig.add_subplot(gs[0, 0])
    ax_oir = fig.add_subplot(gs[1, 0])
    ax_cb = fig.add_subplot(gs[:, 1])

    patches = dataset['patches']
    cmap = perceptual_rainbow_16.mpl_colormap
    radius_normalizer = mpl.colors.Normalize(vmin=0, vmax=22, clip=True)
    r_mapper = mpl.cm.ScalarMappable(norm=radius_normalizer, cmap=cmap)
    radii = []
    for patch_name, patch_group in patches.items():
        for fit_key, ax in zip(('lewis', 'oir_all'), (ax_ms, ax_oir)):
            r_kpc = patch_group.attrs['r_kpc']
            radii.append(r_kpc)
            # logage, sfr = marginalize_metallicity(patch_group, fit_key)
            t = np.array(patch_group['sfh_marginal'][fit_key])
            sfr = t['sfr']
            logage = t['log(age)']
            scaled_sfr = scale_sfr(sfr, patch_group)
            ax.plot(logage, scaled_sfr, '-', lw=0.5,
                    c=r_mapper.to_rgba(r_kpc, alpha=0.5))
    r_mapper.set_array(np.array(radii))
    cbar = fig.colorbar(r_mapper, cax=ax_cb, orientation='vertical')
    cbar.set_label(r'$R_\mathrm{maj}~(\mathrm{kpc})$')
    for ax in (ax_ms, ax_oir):
        ax.set_xlim(6.4, 10.2)
        ax.set_ylabel(SFR_LABEL)
        ax.set_ylim(-5., 10.)
        ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(base=1))
    for tl in ax_ms.get_xmajorticklabels():
        tl.set_visible(False)
    ax_oir.set_xlabel(r'$\log_{10}(A~\mathrm{yr}^{-1})$')
    ax_ms.text(0.1, 0.9, 'ACS-MS', transform=ax_ms.transAxes)
    ax_oir.text(0.1, 0.9, 'OIR-ALL', transform=ax_oir.transAxes)
    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def plot_mean_age_map(dataset, plot_path):
    """Plot maps of the mean stellar age in the patches."""
    fig, canvas, ax_ms, ax_oir, ax_cb = setup_plane_comp_axes()

    ra = dataset['sfh_table']['ra'][:]
    dec = dataset['sfh_table']['dec'][:]

    cmap = perceptual_rainbow_16.mpl_colormap
    normalizer = mpl.colors.Normalize(vmin=0, vmax=8, clip=True)

    for ax, fit_key in zip([ax_ms, ax_oir], ['lewis', 'oir_all']):
        mean_age = dataset['sfh_table']['mean_age_{0}'.format(fit_key)]
        mapper = ax.scatter(ra, dec, c=mean_age, norm=normalizer, cmap=cmap,
                            edgecolors='None', s=50,
                            transform=ax.get_transform('world'))

    cbar = fig.colorbar(mapper, cax=ax_cb, orientation='vertical')
    cbar.set_label(r'$\langle A \rangle$ (Gyr)')
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def plot_single_mean_age_map(dataset, plane_name, plot_path):
    fig, canvas, ax, ax_cb = setup_plane_axes()

    ra = dataset['sfh_table']['ra'][:]
    dec = dataset['sfh_table']['dec'][:]

    cmap = perceptual_rainbow_16.mpl_colormap
    normalizer = mpl.colors.Normalize(vmin=1, vmax=8, clip=True)

    mean_age = dataset['sfh_table']['mean_age_{0}'.format(plane_name)]
    mapper = ax.scatter(ra, dec, c=mean_age, norm=normalizer, cmap=cmap,
                        edgecolors='None', s=45,
                        transform=ax.get_transform('world'))

    cbar = fig.colorbar(mapper, cax=ax_cb, orientation='vertical')
    cbar.set_label(r'$\langle A \rangle$ (Gyr)')
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def plot_epoch_sfr_map_vertical(dataset, fit_key, plot_path):
    ages = [10, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
            2000, 3000, 4000, 5000]
    nx = 4
    ny = 4
    assert len(ages) == nx * ny

    basemap = load_galex_map()

    fig = Figure(figsize=(6.5, 7.5), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(ny, nx + 1,
                           left=0.12, right=0.9, bottom=0.05, top=0.98,
                           wspace=0.05, hspace=0.05,
                           width_ratios=[1] * nx + [0.1],
                           height_ratios=None)
    ax_cb = fig.add_subplot(gs[:, nx])
    axes = {}
    for i, age in enumerate(ages):
        iy = i / nx
        ix = i % nx
        ax = setup_galex_axes(fig, gs[iy, ix], basemap)
        if ix > 0:
            ax.coords[1].ticklabels.set_visible(False)
        if iy < (ny - 1):
            ax.coords[0].ticklabels.set_visible(False)
        axes[age] = ax

    for age, ax in axes.iteritems():
        cmap = perceptual_rainbow_16.mpl_colormap
        normalizer = mpl.colors.Normalize(vmin=0, vmax=5, clip=True)
        ra, dec, sfr_val = get_scaled_sfr_values(dataset, fit_key,
                                                 age * 10. ** 6)
        mapper = ax.scatter(ra, dec, c=sfr_val,
                            norm=normalizer,
                            cmap=cmap,
                            edgecolors='None', s=16,
                            transform=ax.get_transform('world'))
        ax.text(0.95, 0.95, '{0:d} Myr'.format(int(age)),
                ha='right', va='top', transform=ax.transAxes)

    cbar = fig.colorbar(mapper,
                        cax=ax_cb, orientation='vertical')
    cbar.set_label(SFR_LABEL)

    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def plot_major_ax_sfr(dataset, plot_path):
    age_spans = [(0, 25), (25, 50), (50, 79), (79, 100), (100, 158),
                 (158, 200), (200, 251), (251, 316), (316, 398), (0, 400)]
    palette = palettable.cubehelix.Cubehelix.make(
        start_hue=240., end_hue=-300., min_sat=1., max_sat=2.5,
        min_light=0.3, max_light=0.8, gamma=.9, n=len(age_spans) - 1)
    colors = list(palette.mpl_colors) + ['k']
    labels = ['{0} - {1} Myr'.format(*a) for a in age_spans] + ['400 Myr Mean']

    basemap = load_galex_map()
    fig = Figure(figsize=(6.5, 3.0), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(1, 3,
                           left=0.1, right=0.97, bottom=0.15, top=0.95,
                           wspace=0.1, hspace=None,
                           width_ratios=(1, 1, 0.3), height_ratios=None)
    ax_ms = fig.add_subplot(gs[0])
    ax_oir = fig.add_subplot(gs[1])
    for ax in (ax_ms, ax_oir):
        ax.set_xlabel(r'$R_\mathrm{maj}~(\mathrm{kpc})$')
        ax.set_xlim(0, 20.)
        ax.set_ylim(-2, 4)
    ax_ms.set_ylabel(SFR_LABEL)
    for tl in ax_oir.get_ymajorticklabels():
        tl.set_visible(False)
    ax_ms.text(0.9, 0.9, 'ACS-MS', ha='right', va='top',
               transform=ax_ms.transAxes)
    ax_oir.text(0.9, 0.9, 'OIR-ALL', ha='right', va='top',
                transform=ax_oir.transAxes)

    ax_map = setup_galex_axes(fig, gs[2], basemap)
    plot_patch_footprints(ax_map)
    patch_keys = majoraxplot.select_patches(dataset)
    majoraxplot.plot_highlighted_patches(dataset, patch_keys, ax_map)
    ax_map.coords[0].ticklabels.set_visible(False)
    ax_map.coords[1].ticklabels.set_visible(False)

    r_grid, binned_patches = majoraxplot.bin_patches_radially(dataset,
                                                              patch_keys)
    for fit_key, ax in zip(('lewis', 'oir_all'), (ax_ms, ax_oir)):
        for (age_min, age_max), c, label in zip(age_spans, colors, labels):
            sfr = [majoraxplot.compute_sfr_in_span(dataset, patches, fit_key,
                                                   age_min, age_max)
                   for patches in binned_patches]
            ax.plot(r_grid, sfr, c=c, label=label)

    ax_oir.legend(loc='lower left', frameon=False, ncol=2, fontsize=6)

    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def plot_major_ax_sfr_linear(dataset, plot_path):
    age_spans = [(0, 25), (25, 50), (50, 79), (79, 100), (100, 158),
                 (158, 200), (200, 251), (251, 316), (316, 398), (0, 400)]
    palette = palettable.cubehelix.Cubehelix.make(
        start_hue=240., end_hue=-300., min_sat=1., max_sat=2.5,
        min_light=0.3, max_light=0.8, gamma=.9, n=len(age_spans) - 1)
    colors = list(palette.mpl_colors) + ['k']
    labels = ['{0} - {1} Myr'.format(*a) for a in age_spans] + ['400 Myr Mean']

    basemap = load_galex_map()
    fig = Figure(figsize=(6.5, 3.0), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(1, 3,
                           left=0.12, right=0.97, bottom=0.15, top=0.95,
                           wspace=0.1, hspace=None,
                           width_ratios=(1, 1, 0.3), height_ratios=None)
    ax_ms = fig.add_subplot(gs[0])
    ax_oir = fig.add_subplot(gs[1])
    for ax in (ax_ms, ax_oir):
        ax.set_xlabel(r'$R_\mathrm{maj}~(\mathrm{kpc})$')
        ax.set_xlim(0, 20.)
        ax.set_ylim(-5, 4000)
    ax_ms.set_ylabel(LIN_SFR_LABEL)
    for tl in ax_oir.get_ymajorticklabels():
        tl.set_visible(False)
    ax_ms.text(0.9, 0.9, 'ACS-MS', ha='right', va='top',
               transform=ax_ms.transAxes)
    ax_oir.text(0.9, 0.9, 'OIR-ALL', ha='right', va='top',
                transform=ax_oir.transAxes)

    ax_map = setup_galex_axes(fig, gs[2], basemap)
    plot_patch_footprints(ax_map)
    patch_keys = majoraxplot.select_patches(dataset)
    majoraxplot.plot_highlighted_patches(dataset, patch_keys, ax_map)
    ax_map.coords[0].ticklabels.set_visible(False)
    ax_map.coords[1].ticklabels.set_visible(False)

    r_grid, binned_patches = majoraxplot.bin_patches_radially(dataset,
                                                              patch_keys)
    for fit_key, ax in zip(('lewis', 'oir_all'), (ax_ms, ax_oir)):
        for (age_min, age_max), c, label in zip(age_spans, colors, labels):
            sfr = [majoraxplot.compute_sfr_in_span(dataset, patches, fit_key,
                                                   age_min, age_max,
                                                   lin_scale=True)
                   for patches in binned_patches]
            ax.plot(r_grid, sfr, c=c, label=label)

    ax_ms.legend(loc='center left', frameon=False, ncol=2, fontsize=6)

    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def plot_major_ax_cumulative_mass_unnorm(dataset, plot_path):
    data = _prep_cumulative_mass_dataset(dataset)

    fig = Figure(figsize=(6, 4), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(1, 1,
                           left=0.15, right=0.95, bottom=0.15, top=0.95,
                           wspace=None, hspace=None,
                           width_ratios=None, height_ratios=None)
    ax = fig.add_subplot(gs[0])
    ax.set_xlabel(r'$\log(A~\mathrm{yr}^{-1})$')
    ax.set_ylabel(r'$M(t_\mathrm{lookback}<A) [M_\odot]$')
    for patch in data:
        ax.plot(patch['logage'], patch['cmass'], ls='-', lw=1., c='k')

    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def plot_major_ax_cumulative_mass_norm(dataset, plot_path):
    data = _prep_cumulative_mass_dataset(dataset)
    # patches = dataset['patches']
    radii = np.array([p['r_kpc'] for p in data])
    srt = np.argsort(radii)[::-1]  # sort outside in for better visibility
    print srt
    radii = radii[srt]
    data = [data[i] for i in srt]

    fig = Figure(figsize=(6, 4), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(1, 3,
                           left=0.1, right=0.85, bottom=0.15, top=0.95,
                           wspace=0.02, hspace=None,
                           width_ratios=(0.3, 1, 0.05), height_ratios=None)
    ax_log = fig.add_subplot(gs[0])
    ax_lin = fig.add_subplot(gs[1])
    ax_cb = fig.add_subplot(gs[2])

    cmap = perceptual_rainbow_16.mpl_colormap
    radius_normalizer = mpl.colors.Normalize(vmin=0, vmax=22, clip=True)
    r_mapper = mpl.cm.ScalarMappable(norm=radius_normalizer, cmap=cmap)

    ax_log.set_xlabel(r'$\log(A~\mathrm{yr}^{-1})$')
    ax_lin.set_xlabel(r'$A$ (Gyr)')
    ax_log.set_ylabel(r'$M(t_\mathrm{lookback}>A) / \sum M$')
    # No tick labels for linear axes vertical axis
    for tl in ax_lin.get_ymajorticklabels():
        tl.set_visible(False)
    ax_log.xaxis.set_major_locator(mpl.ticker.MultipleLocator(base=1))

    for patch in data:
        normalized_mass = patch['cmass'] / patch['cmass'].max()
        r_kpc = patch['r_kpc']
        plot_args = dict(ls='-',
                         lw=1.,
                         c=r_mapper.to_rgba(r_kpc, alpha=1.))
        ax_log.plot(patch['logage'],
                    normalized_mass,
                    **plot_args)
        ax_lin.plot(10. ** (patch['logage'] - 9.),
                    normalized_mass,
                    **plot_args)

    r_mapper.set_array(np.array(radii))
    cbar = fig.colorbar(r_mapper, cax=ax_cb, orientation='vertical')
    cbar.set_label(r'$R_\mathrm{maj}~(\mathrm{kpc})$')

    ax_log.set_ylim(0., 1.)
    ax_lin.set_ylim(0., 1.)

    ax_log.set_xlim(6.5, 9.)
    ax_lin.set_xlim(1., 13.)

    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def _prep_cumulative_mass_dataset(dataset):
    patch_keys = majoraxplot.select_patches(dataset)
    patches = dataset['patches']
    data = []
    for k in patch_keys:
        patch = patches[k]
        sfh = np.array(patch['sfh_marginal']['oir_all'])
        r_kpc = patch.attrs['r_kpc']
        area = patch.attrs['area_proj'] \
            / np.cos(77.5 * np.pi / 180.) / 1e3 / 1e3  # kpc^2
        cumulative_mass_kpc2 = np.cumsum(sfh['mass'][::-1] / area)[::-1]
        logage = sfh['log(age)']
        data.append({'r_kpc': r_kpc,
                     'cmass': cumulative_mass_kpc2,
                     'logage': logage})
    return data


def plot_mass_accumulation_age_map(dataset, plot_path, fit_key='oir_all'):
    """Similar to plot_epoch_sfr_map_vertical, but plots age when each
    patch reaches a threshold fraction of mass accumulation.
    """
    fractions = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    nx = 4
    ny = 2
    assert len(fractions) == nx * ny

    basemap = load_galex_map()

    fig = Figure(figsize=(6.5, 4.5), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(ny, nx + 1,
                           left=0.12, right=0.9, bottom=0.05, top=0.93,
                           wspace=0.05, hspace=0.1,
                           width_ratios=[1] * nx + [0.1],
                           height_ratios=None)
    ax_cb = fig.add_subplot(gs[:, nx])
    axes = []
    for i, mass_frac in enumerate(fractions):
        iy = i / nx
        ix = i % nx
        ax = setup_galex_axes(fig, gs[iy, ix], basemap)
        if ix > 0:
            ax.coords[1].ticklabels.set_visible(False)
        if iy < (ny - 1):
            ax.coords[0].ticklabels.set_visible(False)
        axes.append(ax)

    for mass_frac, ax in zip(fractions, axes):
        cmap = perceptual_rainbow_16.mpl_colormap
        normalizer = mpl.colors.Normalize(vmin=0, vmax=6, clip=True)
        ra, dec, logage_at_mass_frac = _compute_ages_at_mass_frac(dataset,
                                                                  mass_frac,
                                                                  fit_key)
        age_gyr = 10. ** (logage_at_mass_frac - 9.)
        mapper = ax.scatter(ra, dec, c=age_gyr,
                            norm=normalizer,
                            cmap=cmap,
                            edgecolors='None', s=16,
                            transform=ax.get_transform('world'))
        label = r'$\frac{{M(t_\mathrm{{L}}>A)}}{{\sum M}} ' \
            '\geq {0:.1f}$'.format(mass_frac)
        ax.text(0.5, 1.02, label,
                ha='center', va='bottom', size=10, transform=ax.transAxes,
                backgroundcolor='w')

    cbar = fig.colorbar(mapper,
                        cax=ax_cb, orientation='vertical')
    cbar.set_label(r'$A$ (Gyr)')

    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def _compute_ages_at_mass_frac(dataset, frac_threshold, fit_key):
    ra, dec, thresh_logage = [], [], []

    for k in dataset['patches']:
        p = dataset['patches'][k]
        print k
        print p
        sfh = np.array(p['sfh_marginal'][fit_key])
        ra.append(p.attrs['ra0'])
        dec.append(p.attrs['dec0'])
        area = p.attrs['area_proj'] \
            / np.cos(77.5 * np.pi / 180.) / 1e3 / 1e3  # kpc^2
        cumulative_mass_kpc2 = np.cumsum(sfh['mass'][::-1] / area)[::-1]
        cmass_norm = cumulative_mass_kpc2 / cumulative_mass_kpc2.max()
        logage = sfh['log(age)']
        print frac_threshold, cmass_norm
        try:
            arg_below_fraction = np.argwhere(
                cmass_norm <= frac_threshold)[0][0]
        except IndexError:
            # No age when this mass fraction was recorded;
            # default to oldest age
            arg_below_fraction = len(logage) - 1
        thresh_logage.append(logage[arg_below_fraction])

    return np.array(ra), np.array(dec), np.array(thresh_logage)


if __name__ == '__main__':
    main()
