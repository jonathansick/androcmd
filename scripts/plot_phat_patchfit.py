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


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('hdf5_path')
    parser.add_argument('--radial-sfh-points', default=None)
    parser.add_argument('--radial-sfh-intervals', default=None)
    parser.add_argument('--rchi-hist', default=None)
    parser.add_argument('--sfh-lines', default=None)
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
    ax.set_xlim(0., 30.)
    ax.set_ylim(0., 14.)
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
    radius_normalizer = mpl.colors.Normalize(vmin=0, vmax=25, clip=True)
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
    ax_oir.set_xlabel(r'$\log(A~\mathrm{Gyr}^{-1})$')
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


if __name__ == '__main__':
    main()
