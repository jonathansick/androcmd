#!/usr/bin/env python
# encoding: utf-8
"""
Plotting and analysis of the PHAT patch fitting.

Use the fetch_patch_results.py script to build a patch dataset.

2015-06-17 - Created by Jonathan Sick
"""

import argparse
import h5py

# import matplotlib as mpl
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.gridspec as gridspec


def main():
    args = parse_args()

    dataset = h5py.File(args.hdf5_path, 'r')

    if args.radial_points is not None:
        plot_radial_points(dataset, args.radial_points)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('hdf5_path')
    parser.add_argument('--radial-points', default=None)
    return parser.parse_args()


def plot_radial_points(dataset, plot_path):
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
    print R
    print mean_age
    ax.scatter(R, mean_age, s=3, edgecolors='None', facecolors='dodgerblue')
    ax.scatter(R, mean_age_ms, s=3, edgecolors='None', facecolors='firebrick')
    ax.set_xlim(0., 30.)
    ax.set_ylim(0., 14.)
    ax.set_xlabel(r'$R_\mathrm{maj}$ (kpc)')
    ax.set_ylabel(r'$\langle A \rangle$ (Gyr)')
    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


if __name__ == '__main__':
    main()