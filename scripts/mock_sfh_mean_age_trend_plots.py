#!/usr/bin/env python
# encoding: utf-8
"""
Plot mean age recovery accuracy vs SSP age or tau for different SFH fields.
"""

import os

import numpy as np
from astropy.table import Table
import h5py

import matplotlib as mpl
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.gridspec as gridspec
# from palettable.cubehelix import Cubehelix
# from palettable.colorbrewer.qualitative import Set1_4, Set1_5
from palettable.colorbrewer.qualitative import Dark2_6
# from palettable.colorbrewer.sequential import YlGnBu_4


def main():
    plot_ssp_mean_age_accuracy('real_mock_ssp_mean_age_accuracy')
    plot_ssp_mean_age_accuracy(
        'mock_ssp_mean_age_accuracy',
        experiments=('m3', 'm4', 'm5', 'm6', 'ideal6'),
        labels=('\#3', '\#4', '\#5', '\#6', 'No Errors'),
        colors=Dark2_6.mpl_colors)
    plot_ideal_mean_age_accuracy('mock_ssp_ideal_mean_age_accuracy')
    plot_ssp_mean_age_percent_accuracy(
        'real_mock_ssp_mean_age_percent_accuracy')
    plot_tau_mean_age_accuracy(
        'real_mock_tau_mean_age_accuracy')


def plot_ssp_mean_age_accuracy(plot_path,
                               experiments=('m3', 'm4', 'm5', 'm6'),
                               labels=('\#3', '\#4', '\#5', '\#6'),
                               colors=Dark2_6.mpl_colors):
    fig = Figure(figsize=(6.5, 3.5), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(1, 2,
                           left=0.1, right=0.95, bottom=0.15, top=0.95,
                           wspace=0.1, hspace=None,
                           width_ratios=None, height_ratios=None)
    ax_lewis = fig.add_subplot(gs[0])
    ax_oirall = fig.add_subplot(gs[1])
    for plane, ax in zip(['lewis', 'oir_all'], [ax_lewis, ax_oirall]):
        for experiment, label, color in zip(experiments, labels, colors[:4]):
            t = ssp_mean_age_table(experiment=experiment)
            age_diff = t['{0}_mean_age'.format(plane)] - t['mock_age']
            age_sigma = t['{0}_mean_age_sigma'.format(plane)]
            ax.plot(t['mock_age'],
                    age_diff,
                    c=color, label=label, lw=1.5)
            ax.plot(t['mock_age'],
                    np.array(-age_sigma),
                    c=color, lw=0.5, label=None, ls='--')
            ax.plot(t['mock_age'],
                    np.array(age_sigma),
                    c=color, lw=0.5, label=None, ls='--')
        ax.set_xlabel(r'$\langle A \rangle_\mathrm{mock}$ (Gyr)')
        ax.set_ylim(-10, 10)
    ax_lewis.set_ylabel(
        r'$\langle A \rangle_\mathrm{fit} - \langle A \rangle_\mathrm{mock}$ (Gyr)')  # NOQA
    for tl in ax_oirall.get_ymajorticklabels():
        tl.set_visible(False)
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(base=1))
    ax_lewis.text(0.1, 0.9, 'ACS-MS', transform=ax_lewis.transAxes)
    ax_oirall.text(0.1, 0.9, 'OIR-ALL', transform=ax_oirall.transAxes)
    ax_lewis.set_xlim(0., 2.)
    ax_lewis.legend(loc='lower left')
    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def plot_ideal_mean_age_accuracy(plot_path):
    colors = Dark2_6.mpl_colors
    print colors

    real_experiments = ['m3', 'm4', 'm5', 'm6']
    real_labels = [r'\#3', r'\#4', r'\#5', r'\#6']
    real_colors = colors[:4]

    ideal_experiment = 'idealall'
    ideal_planes = ['oir_all', 'oir_all_28']
    ideal_labels = ['Errorless OIR-ALL',
                    'Errorless OIR-ALL-28']
    ideal_colors = colors[4:]
    print ideal_colors

    fig = Figure(figsize=(6.5, 3.5), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(1, 1,
                           left=0.1, right=0.95, bottom=0.15, top=0.95,
                           wspace=0.1, hspace=None,
                           width_ratios=None, height_ratios=None)
    ax = fig.add_subplot(gs[0])

    itr = zip(real_experiments, real_labels, real_colors)
    for experiment, label, color in itr:
        t = ssp_mean_age_table(experiment=experiment)
        age_diff = t['{0}_mean_age'.format('oir_all')] - t['mock_age']
        # age_sigma = t['{0}_mean_age_sigma'.format('oir_all')]
        print experiment, color
        ax.plot(t['mock_age'],
                age_diff,
                c=color, label=label, lw=3)
                # dashes=(5, 7.5),
                # dash_joinstyle='round', dash_capstyle='round')
        # ax.plot(t['mock_age'],
        #         -age_sigma,
        #         c=color, lw=0.5, ls='--')
        # ax.plot(t['mock_age'],
        #         age_sigma,
        #         c=color, lw=0.5, ls='--')

    for plane, label, color in zip(ideal_planes, ideal_labels, ideal_colors):
        t = ssp_mean_age_table(experiment=ideal_experiment,
                               plane_keys=ideal_planes)
        age_diff = t['{0}_mean_age'.format(plane)] - t['mock_age']
        # age_sigma = t['{0}_mean_age_sigma'.format(plane)]
        print label, color
        ax.plot(t['mock_age'],
                age_diff,
                c=color, label=label, lw=3,
                dashes=(5, 7.5),
                dash_joinstyle='round', dash_capstyle='round')

    ax.set_xlabel(r'$\langle A \rangle_\mathrm{mock}$ (Gyr)')
    ax.set_ylim(-10., 10.)

    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(base=1))

    ax.set_ylabel(
        r'$\langle A \rangle_\mathrm{fit} - \langle A \rangle_\mathrm{mock}$ (Gyr)')  # NOQA
    ax.legend(loc='lower left', fontsize=7, handlelength=5)
    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def plot_ssp_mean_age_percent_accuracy(plot_path):
    experiments = ['m3', 'm4', 'm5', 'm6']
    labels = ['\#3', '\#4', '\#5', '\#6']
    colors = Dark2_6.mpl_colors

    fig = Figure(figsize=(6.5, 3.5), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(1, 2,
                           left=0.12, right=0.95, bottom=0.15, top=0.95,
                           wspace=0.1, hspace=None,
                           width_ratios=None, height_ratios=None)
    ax_lewis = fig.add_subplot(gs[0])
    ax_oirall = fig.add_subplot(gs[1])
    for plane, ax in zip(['lewis', 'oir_all'], [ax_lewis, ax_oirall]):
        for experiment, label, color in zip(experiments, labels, colors):
            t = ssp_mean_age_table(experiment=experiment)
            age_diff = (t['{0}_mean_age'.format(plane)] - t['mock_age']) \
                / t['mock_age'] * 100.
            age_sigma = t['{0}_mean_age_sigma'.format(plane)] \
                / t['mock_age'] * 100.
            ax.plot(t['mock_age'],
                    age_diff,
                    c=color, label=label, lw=1.5)
            ax.plot(t['mock_age'],
                    -age_sigma,
                    c=color, lw=0.5, ls='--')
            ax.plot(t['mock_age'],
                    age_sigma,
                    c=color, lw=0.5, ls='--')
        ax.set_xlabel(r'$\langle A \rangle_\mathrm{mock}$ (Gyr)')
        ax.set_ylim(-400., 400.)
    ax_lewis.set_ylabel(
        r'$\langle A \rangle_\mathrm{fit} - \langle A \rangle_\mathrm{mock}$ (\%)')  # NOQA
    for tl in ax_oirall.get_ymajorticklabels():
        tl.set_visible(False)
    ax_lewis.text(0.1, 0.9, 'ACS-MS', transform=ax_lewis.transAxes)
    ax_oirall.text(0.1, 0.9, 'OIR-ALL', transform=ax_oirall.transAxes)
    ax_lewis.legend(loc='upper right')
    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def plot_tau_mean_age_accuracy(plot_path):
    experiments = ['m3', 'm4', 'm5', 'm6']
    labels = ['\#3', '\#4', '\#5', '\#6']
    colors = Dark2_6.mpl_colors

    fig = Figure(figsize=(6.5, 3.5), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(1, 2,
                           left=0.1, right=0.95, bottom=0.15, top=0.95,
                           wspace=0.1, hspace=None,
                           width_ratios=None, height_ratios=None)
    ax_lewis = fig.add_subplot(gs[0])
    ax_oirall = fig.add_subplot(gs[1])
    for plane, ax in zip(['lewis', 'oir_all'], [ax_lewis, ax_oirall]):
        for experiment, label, color in zip(experiments, labels, colors):
            t = tau_mean_age_table(experiment=experiment)
            age_diff = t['{0}_mean_age'.format(plane)] - t['mock_age']
            age_sigma = t['{0}_mean_age_sigma'.format(plane)]
            ax.plot(np.log10(t['mock_tau']),
                    age_diff,
                    c=color, label=label, lw=1.5)
            ax.plot(np.log10(t['mock_age']),
                    -age_sigma,
                    c=color, lw=0.5, ls='--')
            ax.plot(np.log10(t['mock_age']),
                    age_sigma,
                    c=color, lw=0.5, ls='--')
        ax.set_xlabel(r'$\log_{10} \tau_\mathrm{mock}$')
        ax.set_ylim(-10, 10)
    ax_lewis.set_ylabel(
        r'$\langle A \rangle_\mathrm{fit} - \langle A \rangle_\mathrm{mock}$ (Gyr)')  # NOQA
    for tl in ax_oirall.get_ymajorticklabels():
        tl.set_visible(False)
    ax_lewis.text(0.1, 0.9, 'ACS-MS', transform=ax_lewis.transAxes)
    ax_oirall.text(0.1, 0.9, 'OIR-ALL', transform=ax_oirall.transAxes)
    ax_lewis.legend(loc='lower right')
    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def ssp_mean_age_table(experiment='m3', plane_keys=('lewis', 'oir_all')):
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
    cols = [mock_mean_ages] + [fit_mean_ages[p] for p in plane_keys] \
        + [fit_mean_age_sigmas[p] for p in plane_keys]
    colnames = ['mock_age'] + ['{}_mean_age'.format(p) for p in plane_keys] \
        + ['{}_mean_age_sigma'.format(p) for p in plane_keys]
    return Table(cols, names=colnames)


def tau_mean_age_table(experiment='m3'):
    """Make a table of tau age (Myr), recovered age (Myr)."""
    dirname = os.path.join(os.getenv('STARFISH'), experiment)
    h5path = os.path.join(dirname, experiment + '.hdf5')
    hdf5 = h5py.File(h5path, 'r')
    dataset = hdf5['mocksfh']

    mock_list = ['tau_0.1_solar',
                 'tau_0.5_solar',
                 'tau_1.0_solar',
                 'tau_5.0_solar',
                 'tau_10.0_solar',
                 'tau_20.0_solar',
                 'tau_50.0_solar',
                 'tau_100.0_solar']
    mock_tau = np.array([0.1, 0.5, 1.0, 5.0, 10.0, 20.0, 50.0, 100.0])

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
    cols = [mock_mean_ages, mock_tau,
            fit_mean_ages['lewis'], fit_mean_age_sigmas['lewis'],
            fit_mean_ages['oir_all'], fit_mean_age_sigmas['oir_all']]
    colnames = ['mock_age', 'mock_tau',
                'lewis_mean_age', 'lewis_mean_age_sigma',
                'oir_all_mean_age', 'oir_all_mean_age_sigma']
    return Table(cols, names=colnames)


if __name__ == '__main__':
    main()
