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
import wcsaxes
import astropy.io.fits

from androcmd.phatpatchfit import load_field_footprints


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


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('hdf5_path')
    parser.add_argument('--radial-sfh-points', default=None)
    parser.add_argument('--radial-sfh-intervals', default=None)
    parser.add_argument('--rchi-hist', default=None)
    parser.add_argument('--sfh-lines', default=None)
    parser.add_argument('--mean-sfr-map', default=None)
    parser.add_argument('--mean-age-map', default=None)
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


def marginalize_metallicity(sfh_table):
    t = np.empty(len(sfh_table), dtype=sfh_table.dtype)
    sfh_table.read_direct(t, source_sel=None, dest_sel=None)
    age_vals = np.unique(t['log(age)'])
    s = np.argsort(age_vals)
    age_vals = age_vals[s]
    A = []
    sfr = []
    for i, age_val in enumerate(age_vals):
        tt = t[t['log(age)'] == age_val]
        bin_sfr = np.sum(tt['sfr'])
        A.append(age_val)
        sfr.append(bin_sfr)
    srt = np.argsort(A)
    A = np.array(A)
    sfr = np.array(sfr)
    A = A[srt]
    sfr = sfr[srt]
    return A, sfr


def plot_mean_sfr_map(dataset, plot_path):
    """Plot maps of the mean star formation rate in the patches."""
    fig, canvas, ax_ms, ax_oir, ax_cb = create_wcs_axes_galex()

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
    fig, canvas, ax_ms, ax_oir, ax_cb = create_wcs_axes_galex()

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


def create_wcs_axes(ref_path='m31_80.fits'):
    fig = Figure(figsize=(6, 3.0), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(1, 3,
                           left=0.08, right=0.85, bottom=0.15, top=0.95,
                           wspace=0.05, hspace=None,
                           width_ratios=(1, 1, 0.08), height_ratios=None)

    with astropy.io.fits.open(ref_path) as f:
        header = f[0].header
        base_image = f[0].data
        wcs = wcsaxes.WCS(header)

    ax_ms = fig.add_subplot(gs[0], projection=wcs)
    ax_oir = fig.add_subplot(gs[1], projection=wcs)
    ax_cb = fig.add_subplot(gs[2])

    for ax in (ax_oir, ax_ms):
        # ax.invert_xaxis()
        ax.set_xlim(-0.5, base_image.shape[1] - 0.5)
        ax.set_ylim(-0.5, base_image.shape[0] - 0.5)
        ax.imshow(np.log10(base_image),
                  cmap=mpl.cm.gray_r, vmin=0.3, vmax=0.6,
                  zorder=-10,
                  origin='lower')
        ax.set_xlim(2000, 9500)
        ax.coords[1].set_major_formatter('d.d')
        ax.coords[0].set_major_formatter('hh:mm')
    # ax_oir.coords[1].ticklabels.set_visible(False)
    ax_oir.coords[1].set_ticklabel_position('')
    ax_ms.text(0.1, 0.9, 'ACS-MS', transform=ax_ms.transAxes, ha='left',
               zorder=10)
    ax_oir.text(0.1, 0.9, 'OIR-ALL', transform=ax_oir.transAxes, ha='left',
                zorder=10)

    return fig, canvas, ax_ms, ax_oir, ax_cb


def create_wcs_axes_galex(ref_path='h_m31-nd-int.fits'):
    fig = Figure(figsize=(6, 3.5), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(1, 3,
                           left=0.08, right=0.85, bottom=0.1, top=0.90,
                           wspace=0.05, hspace=None,
                           width_ratios=(1, 1, 0.08), height_ratios=None)

    with astropy.io.fits.open(ref_path) as f:
        header = f[0].header
        base_image = f[0].data
        wcs = wcsaxes.WCS(header)

    ax_ms = fig.add_subplot(gs[0], projection=wcs)
    ax_oir = fig.add_subplot(gs[1], projection=wcs)
    ax_cb = fig.add_subplot(gs[2])

    for ax in (ax_oir, ax_ms):
        # ax.invert_xaxis()
        ax.set_xlim(-0.5, base_image.shape[1] - 0.5)
        ax.set_ylim(-0.5, base_image.shape[0] - 0.5)
        ax.imshow(np.log10(base_image),
                  cmap=mpl.cm.gray_r, vmin=-2., vmax=-1,
                  zorder=-10,
                  origin='lower')
        ax.set_xlim(500, 3000)
        ax.set_ylim(2800, 6189)
        ax.coords[1].set_major_formatter('d.d')
        ax.coords[0].set_major_formatter('hh:mm')
        ax.coords[0].set_separator(('h', "'", '"'))

        # Plot phat footprints
        for footprint in load_field_footprints():
            patch = mpl.patches.Polygon(footprint, closed=True,
                                        transform=ax.get_transform('world'),
                                        facecolor='None', alpha=0.2,
                                        edgecolor='k', lw=0.5)
            ax.add_patch(patch)

    ax_oir.coords[1].ticklabels.set_visible(False)
    ax_ms.text(0.1, 0.1, 'ACS-MS', transform=ax_ms.transAxes, ha='left',
               zorder=10)
    ax_oir.text(0.1, 0.1, 'OIR-ALL', transform=ax_oir.transAxes, ha='left',
                zorder=10)

    return fig, canvas, ax_ms, ax_oir, ax_cb


if __name__ == '__main__':
    main()
