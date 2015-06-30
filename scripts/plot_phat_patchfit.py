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
from palettable.cubehelix import perceptual_rainbow_16

from androcmd.phatpatchfit import (load_galex_map,
                                   setup_galex_axes, setup_plane_comp_axes,
                                   marginalize_metallicity,
                                   get_scaled_sfr_values)


def main():
    args = parse_args()

    dataset = h5py.File(args.hdf5_path, 'r')

    if args.radial_sfh_points is not None:
        plot_radial_sfh_points(dataset, args.radial_sfh_points)

    if args.radial_sfh_intervals is not None:
        plot_radial_sfh_intervals(dataset, args.radial_sfh_intervals)

    if args.rchi_hist is not None:
        plot_reduced_chi_hist(dataset, args.rchi_hist)

    if args.sfh_lines is not None:
        plot_sfh_lines(dataset, args.sfh_lines)
        plot_sfh_lines_phi(dataset, args.sfh_lines + "_phi")

    if args.mean_sfr_map is not None:
        plot_mean_sfr_map(dataset, args.mean_sfr_map)

    if args.mean_age_map is not None:
        plot_mean_age_map(dataset, args.mean_age_map)

    if args.epoch_sfr_comp_maps is not None:
        plot_epoch_sfr_comp_maps(dataset, args.epoch_sfr_comp_maps)

    if args.epoch_sfr_maps is not None:
        for fit_key in ['lewis', 'oir_all']:
            plot_epoch_sfr_map_vertical(
                dataset, fit_key,
                '_'.join((args.epoch_sfr_maps, fit_key, 'vert')))


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('hdf5_path')
    parser.add_argument('--radial-sfh-points', default=None)
    parser.add_argument('--radial-sfh-intervals', default=None)
    parser.add_argument('--rchi-hist', default=None)
    parser.add_argument('--sfh-lines', default=None)
    parser.add_argument('--mean-sfr-map', default=None)
    parser.add_argument('--mean-age-map', default=None)
    parser.add_argument(
        '--epoch-sfr-comp-maps', default=None,
        help='Map comparing SFR at epochs for MS and ALL fits')
    parser.add_argument(
        '--epoch-sfr-maps', default=None,
        help='Map comparing SFR at for MS and ALL, individually')
    return parser.parse_args()


def plot_radial_sfh_points(dataset, plot_path):
    fig = Figure(figsize=(3.5, 3.5), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(1, 1,
                           left=0.15, right=0.95, bottom=0.15, top=0.95,
                           wspace=None, hspace=None,
                           width_ratios=None, height_ratios=None)
    ax = fig.add_subplot(gs[0])
    R = dataset['sfh_table']['r_kpc'][:]
    mean_age = dataset['sfh_table']['mean_age_oir_all'][:]
    mean_age_ms = dataset['sfh_table']['mean_age_lewis'][:]
    ax.scatter(R, mean_age, s=3, edgecolors='None', facecolors='firebrick',
               label='OIR-ALL')
    ax.scatter(R, mean_age_ms, s=3, edgecolors='None', facecolors='dodgerblue',
               label='ACS-MS')
    ax.set_xlim(0., 25.)
    ax.set_ylim(0., 12.)
    ax.set_xlabel(r'$R_\mathrm{maj}$ (kpc)')
    ax.set_ylabel(r'$\langle A \rangle$ (Gyr)')
    ax.legend(frameon=True, markerscale=2, ncol=2)
    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def plot_radial_sfh_intervals(dataset, plot_path):
    fig = Figure(figsize=(3.5, 3.5), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(1, 1,
                           left=0.15, right=0.95, bottom=0.15, top=0.95,
                           wspace=None, hspace=None,
                           width_ratios=None, height_ratios=None)
    ax = fig.add_subplot(gs[0])

    R = dataset['sfh_table']['r_kpc'][:]
    age_mean = dataset['sfh_table']['mean_age_oir_all'][:]
    age_25 = dataset['sfh_table']['age_25_oir_all'][:]
    age_75 = dataset['sfh_table']['age_75_oir_all'][:]

    ax.scatter(R, age_mean, s=3, edgecolors='None', facecolors='dodgerblue')
    ax.scatter(R, age_75, s=3, marker='v',
               edgecolors='None', facecolors='dodgerblue')
    ax.scatter(R, age_25, s=3, marker='^',
               edgecolors='None', facecolors='dodgerblue')

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
                           left=0.1, right=0.9, bottom=0.1, top=0.95,
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
            sfh_table = patch_group['sfh'][fit_key]
            logage, sfr = marginalize_metallicity(sfh_table)
            ax.plot(logage, np.log10(sfr), '-', lw=0.5,
                    c=r_mapper.to_rgba(r_kpc, alpha=0.5))
    r_mapper.set_array(np.array(radii))
    cbar = fig.colorbar(r_mapper, cax=ax_cb, orientation='vertical')
    cbar.set_label(r'$R_\mathrm{maj}$')
    for ax in (ax_ms, ax_oir):
        ax.set_xlim(6.4, 10.2)
        ax.set_ylabel(r'$\log(\mathrm{SFR})$')
        ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(base=1))
    for tl in ax_ms.get_xmajorticklabels():
        tl.set_visible(False)
    ax_oir.set_xlabel(r'$\log(A~\mathrm{yr}^{-1})$')
    ax_ms.text(0.1, 0.9, 'ACS-MS', transform=ax_ms.transAxes)
    ax_oir.text(0.1, 0.9, 'OIR-ALL', transform=ax_oir.transAxes)
    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def plot_sfh_lines_phi(dataset, plot_path):
    """Plot SFH for each patch, keyed by PA"""
    fig = Figure(figsize=(6, 6), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(2, 2,
                           left=0.1, right=0.9, bottom=0.1, top=0.95,
                           wspace=0.1, hspace=0.1,
                           width_ratios=(1, 0.1), height_ratios=None)
    ax_ms = fig.add_subplot(gs[0, 0])
    ax_oir = fig.add_subplot(gs[1, 0])
    ax_cb = fig.add_subplot(gs[:, 1])

    patches = dataset['patches']
    cmap = perceptual_rainbow_16.mpl_colormap
    normalizer = mpl.colors.Normalize(vmin=0, vmax=360., clip=True)
    mapper = mpl.cm.ScalarMappable(norm=normalizer, cmap=cmap)
    map_values = []
    for patch_name, patch_group in patches.items():
        for fit_key, ax in zip(('lewis', 'oir_all'), (ax_ms, ax_oir)):
            phi = patch_group.attrs['phi']
            map_values.append(phi)
            sfh_table = patch_group['sfh'][fit_key]
            logage, sfr = marginalize_metallicity(sfh_table)
            ax.plot(logage, np.log10(sfr), '-', lw=0.5,
                    c=mapper.to_rgba(phi, alpha=0.5))
    print min(map_values), max(map_values)
    mapper.set_array(np.array(map_values))
    cbar = fig.colorbar(mapper, cax=ax_cb, orientation='vertical')
    cbar.set_label(r'$\phi$ (degree)')
    for ax in (ax_ms, ax_oir):
        ax.set_xlim(6.4, 10.2)
        ax.set_ylabel(r'$\log(\mathrm{SFR})$')
        ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(base=1))
    for tl in ax_ms.get_xmajorticklabels():
        tl.set_visible(False)
    ax_ms.text(0.1, 0.9, 'ACS-MS', transform=ax_ms.transAxes)
    ax_oir.text(0.1, 0.9, 'OIR-ALL', transform=ax_oir.transAxes)
    ax_oir.set_xlabel(r'$\log(A~\mathrm{yr}^{-1})$')
    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def plot_mean_sfr_map(dataset, plot_path):
    """Plot maps of the mean star formation rate in the patches."""
    fig, canvas, ax_ms, ax_oir, ax_cb = setup_plane_comp_axes()

    cmap = perceptual_rainbow_16.mpl_colormap
    normalizer = mpl.colors.Normalize(vmin=-12, vmax=-4, clip=True)

    patches = dataset['patches']
    for ax, fit_key in zip([ax_ms, ax_oir], ['lewis', 'oir_all']):
        ra = []
        dec = []
        mean_logsfr = []
        for patch_name, patch_group in patches.items():
            sfh_table = patch_group['sfh'][fit_key]
            t = np.empty(len(sfh_table), dtype=sfh_table.dtype)
            sfh_table.read_direct(t, source_sel=None, dest_sel=None)
            mean_logsfr.append(np.log10(t['sfr']).mean())
            ra.append(patch_group.attrs['ra0'])
            dec.append(patch_group.attrs['dec0'])
        mapper = ax.scatter(ra, dec, c=mean_logsfr, norm=normalizer, cmap=cmap,
                            edgecolors='None', s=16,
                            transform=ax.get_transform('world'))

    cbar = fig.colorbar(mapper, cax=ax_cb, orientation='vertical')
    cbar.set_label(r'$\langle \mathrm{SFR} \rangle$')

    canvas.print_figure(plot_path + ".pdf", format="pdf")


def plot_mean_age_map(dataset, plot_path):
    """Plot maps of the mean stellar age in the patches."""
    fig, canvas, ax_ms, ax_oir, ax_cb = setup_plane_comp_axes()

    ra = dataset['sfh_table']['ra'][:]
    dec = dataset['sfh_table']['dec'][:]

    cmap = perceptual_rainbow_16.mpl_colormap
    normalizer = mpl.colors.Normalize(vmin=0, vmax=10, clip=True)

    for ax, fit_key in zip([ax_ms, ax_oir], ['lewis', 'oir_all']):
        mean_age = dataset['sfh_table']['mean_age_{0}'.format(fit_key)]
        mapper = ax.scatter(ra, dec, c=mean_age, norm=normalizer, cmap=cmap,
                            edgecolors='None', s=16,
                            transform=ax.get_transform('world'))

    cbar = fig.colorbar(mapper, cax=ax_cb, orientation='vertical')
    cbar.set_label(r'$\langle A \rangle$ (Gyr)')
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def plot_epoch_sfr_comp_maps(dataset, plot_path):
    fit_keys = ['lewis', 'oir_all']
    ages = [100, 250, 500, 1000]  # Myr

    basemap = load_galex_map()

    fig = Figure(figsize=(9, 5), frameon=False)
    canvas = FigureCanvas(fig)
    nfig = len(ages)
    gs = gridspec.GridSpec(2, nfig + 1,
                           left=0.05, right=0.9, bottom=0.15, top=0.95,
                           wspace=0.05, hspace=0.05,
                           width_ratios=(1, 1, 1, 1, 0.1), height_ratios=None)
    axes_ms = [setup_galex_axes(fig, gs[0, i], basemap)
               for i in range(nfig)]
    axes_oir = [setup_galex_axes(fig, gs[1, i], basemap)
                for i in range(nfig)]
    ax_cb_oir = fig.add_subplot(gs[1, nfig])
    ax_cb_ms = fig.add_subplot(gs[0, nfig])

    for i, (age, ax_ms, ax_oir) in enumerate(zip(ages, axes_ms, axes_oir)):
        for ax in (ax_ms, ax_oir):
            if i > 0:
                ax.coords[1].ticklabels.set_visible(False)
        ax_ms.coords[0].ticklabels.set_visible(False)

        ax_ms.text(0.0, 1.02, '{0:d} Myr'.format(int(age)),
                   ha='left', va='baseline', transform=ax_ms.transAxes)

        # Compute the SFR at age for both ACS-MS and OIR-ALL fits
        mappers = {'lewis': None, 'oir_all': None}
        for fit_key, ax in zip(fit_keys, (ax_ms, ax_oir)):
            cmap = perceptual_rainbow_16.mpl_colormap
            normalizer = mpl.colors.Normalize(vmin=0, vmax=5, clip=True)
            ra, dec, sfr_val = get_scaled_sfr_values(dataset, fit_key,
                                                     age * 10. ** 6)
            mapper = ax.scatter(ra, dec, c=sfr_val,
                                norm=normalizer,
                                cmap=cmap,
                                edgecolors='None', s=16,
                                transform=ax.get_transform('world'))
            mappers[fit_key] = mapper

    sfr_label = r'$\langle \mathrm{SFR} \rangle$ ($\mathrm{M}_\odot' \
                r'~10^{-3}~\mathrm{yr}^{-1}~\mathrm{kpc}^{-2}$)'
    cbar = fig.colorbar(mappers['oir_all'],
                        cax=ax_cb_oir, orientation='vertical')
    cbar.set_label(sfr_label)
    cbar = fig.colorbar(mappers['lewis'], cax=ax_cb_ms,
                        orientation='vertical')
    cbar.set_label(sfr_label)

    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def plot_epoch_sfr_map_vertical(dataset, fit_key, plot_path):
    ages = [10, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    nx = 4
    ny = 3
    assert len(ages) == nx * ny

    basemap = load_galex_map()

    fig = Figure(figsize=(6.5, 6), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(ny, nx + 1,
                           left=0.07, right=0.9, bottom=0.05, top=0.95,
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

    sfr_label = r'$\langle \mathrm{SFR} \rangle$ ($\mathrm{M}_\odot' \
                r'~10^{-3}~\mathrm{yr}^{-1}~\mathrm{kpc}^{-2}$)'
    cbar = fig.colorbar(mapper,
                        cax=ax_cb, orientation='vertical')
    cbar.set_label(sfr_label)

    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


if __name__ == '__main__':
    main()
