#!/usr/bin/env python
# encoding: utf-8
"""
Plot mean age recovery accuracy vs SSP age or tau for different SFH fields.
"""

import os

import numpy as np
from astropy.table import Table
import h5py

# import matplotlib as mpl
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.gridspec as gridspec
# from palettable.cubehelix import Cubehelix
from palettable.colorbrewer.qualitative import Set1_4


def main():
    plot_ssp_mean_age_accuracy('real_mock_ssp_mean_age_accuracy')


def plot_ssp_mean_age_accuracy(plot_path):
    experiments = ['m3', 'm4', 'm5', 'm6']
    labels = ['\#3', '\#4', '\#5', '\#6']
    # palette = Cubehelix.make(start=0.3, rotation=-0.5, n=4)
    # colors = palette.mpl_colors
    colors = Set1_4.mpl_colors

    fig = Figure(figsize=(3.5, 3.5), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(1, 1,
                           left=0.2, right=0.95, bottom=0.15, top=0.95,
                           wspace=None, hspace=None,
                           width_ratios=None, height_ratios=None)
    ax = fig.add_subplot(gs[0])
    for experiment, label, color in zip(experiments, labels, colors):
        t = ssp_mean_age_table(experiment=experiment)
        ax.plot(t['mock_age'], t['oir_all_mean_age'] - t['mock_age'],
                c=color, label=label)
    ax.set_ylabel(
        r'$\langle A \rangle_\mathrm{fit} - \langle A \rangle_\mathrm{mock}$')
    ax.set_xlabel(r'$\langle A \rangle_\mathrm{mock}$')
    ax.legend(loc='lower left')
    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def ssp_mean_age_table(experiment='m3'):
    """Make a table of SSP age (Myr), recovered age (Myr)."""
    dirname = os.path.join(os.getenv('STARFISH'), experiment)
    h5path = os.path.join(dirname, experiment + '.hdf5')
    hdf5 = h5py.File(h5path, 'r')
    dataset = hdf5['mocksfh']

    mock_list = ['ssp_53myr_solar',
                 'ssp_100myr_solar',
                 'ssp_186myr_solar',
                 'ssp_346myr_solar',
                 'ssp_645myr_solar',
                 'ssp_1380myr_solar',
                 'ssp_3467myr_solar',
                 'ssp_5370myr_solar',
                 'ssp_9549myr_solar']

    plane_keys = ['lewis', 'oir_all']
    fit_mean_ages = {k: np.empty(len(mock_list)) for k in plane_keys}
    fit_mean_age_sigmas = {k: np.empty(len(mock_list)) for k in plane_keys}
    mock_mean_ages = np.empty(len(mock_list))
    for i, sfh_key in enumerate(mock_list):
        mock_mean_ages[i] \
            = dataset[sfh_key]['mock_sfh_marginal'].attrs['mean_age'][0]
        for plane_key in plane_keys:
            fit_mean_ages[plane_key][i] \
                = dataset[sfh_key]['sfh'][plane_key].attrs['mean_age'][0]
            fit_mean_age_sigmas[plane_key][i] \
                = dataset[sfh_key]['sfh'][plane_key].attrs['mean_age'][1]
    cols = [mock_mean_ages,
            fit_mean_ages['lewis'], fit_mean_age_sigmas['lewis'],
            fit_mean_ages['oir_all'], fit_mean_age_sigmas['oir_all']]
    colnames = ['mock_age', 'lewis_mean_age', 'lewis_mean_age_sigma',
                'oir_all_mean_age', 'oir_all_mean_age_sigma']
    return Table(cols, names=colnames)


if __name__ == '__main__':
    main()
