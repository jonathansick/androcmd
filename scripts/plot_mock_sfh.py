#!/usr/bin/env python
# encoding: utf-8
"""
Plot mock tests made by the mock_sfh_test.py script.

2015-07-09 - Created by Jonathan Sick
"""

import os
import argparse

import numpy as np
import h5py

import matplotlib as mpl
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.gridspec as gridspec

from palettable.cubehelix import perceptual_rainbow_16
from palettable.colorbrewer.diverging import RdBu_11

from starfisher.sfhplot import plot_single_sfh_line


def main():
    args = parse_args()
    dirname = os.path.join(os.getenv('STARFISH'), args.name)
    h5path = os.path.join(dirname, args.name + '.hdf5')
    hdf5 = h5py.File(h5path, 'r')

    plot_hess_planes(hdf5['mocksfh'], dirname, sfh_list=args.sfh)
    plot_sfhs(hdf5['mocksfh'], dirname, sfh_list=args.sfh)

    hdf5.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('name',
                        help='Test name')
    parser.add_argument('--sfh',
                        help='optional name of sfh to plot',
                        default=None,
                        nargs='*')
    return parser.parse_args()


def plot_sfhs(dataset, base_dir, sfh_list=None):
    if sfh_list is None:
        sfh_list = dataset.keys()
    for sfh_run in sfh_list:
        plot_path = os.path.join(base_dir,
                                 '{0}_sfh'.format(sfh_run))
        plot_sfh(dataset[sfh_run], dataset[sfh_run]['mock_sfh_marginal'],
                 plot_path)


def plot_hess_planes(dataset, base_dir, sfh_list=None):
    if sfh_list is None:
        sfh_list = dataset.keys()
    for sfh_run in sfh_list:
        plane_keys = dataset[sfh_run]['obs_hess'].keys()
        for plane_key in plane_keys:
            plot_path = os.path.join(base_dir,
                                     '{0}_{1}_hess'.format(sfh_run, plane_key))
            _plot_hess(dataset[sfh_run], plane_key, plot_path)


def _plot_hess(dataset, plane_key, plot_path):
    fig = Figure(figsize=(6.5, 3.5), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(2, 4,
                           left=0.08, right=0.97, bottom=0.15, top=0.90,
                           wspace=0.15, hspace=0.,
                           width_ratios=None, height_ratios=(0.1, 1))

    print "dataset.keys()", dataset.keys()
    extent = dataset['fit_hess'][plane_key].attrs['extent']
    origin = dataset['fit_hess'][plane_key].attrs['origin']
    y_label = dataset['fit_hess'][plane_key].attrs['y_label']
    x_label = dataset['fit_hess'][plane_key].attrs['x_label']
    x_span = dataset['fit_hess'][plane_key].attrs['x_span']
    y_span = dataset['fit_hess'][plane_key].attrs['y_span']

    cube_map = perceptual_rainbow_16.mpl_colormap
    div_map = RdBu_11.mpl_colormap

    ax_obs = fig.add_subplot(gs[1, 0])
    ax_model = fig.add_subplot(gs[1, 1])
    ax_chi = fig.add_subplot(gs[1, 2])
    ax_diff = fig.add_subplot(gs[1, 3])
    ax_obs_cb = fig.add_subplot(gs[0, 0])
    ax_model_cb = fig.add_subplot(gs[0, 1])
    ax_chi_cb = fig.add_subplot(gs[0, 2])
    ax_diff_cb = fig.add_subplot(gs[0, 3])

    _imshow = dict(cmap=cube_map,
                   norm=None,
                   aspect='auto',
                   interpolation='none',
                   extent=extent,
                   origin=origin,
                   alpha=None)
    _imshow_diff = dict(_imshow)
    _imshow_diff['cmap'] = div_map

    label_bbox = {'facecolor': 'w',
                  'edgecolor': 'None',
                  'alpha': 0.8,
                  'pad': 5}

    # observations
    hess = np.log10(dataset['obs_hess'][plane_key])
    hess = np.ma.masked_invalid(hess, copy=True)
    im = ax_obs.imshow(hess, vmin=0, vmax=3, **_imshow)
    ax_obs.set_xlabel(x_label)
    ax_obs.set_ylabel(y_label)
    cb = fig.colorbar(im, cax=ax_obs_cb, orientation='horizontal')
    cb.set_label(r"$\log(N_*)$ Mock", bbox=label_bbox)
    cb.ax.xaxis.set_ticks_position('top')
    cb.locator = mpl.ticker.MultipleLocator(1.)
    cb.update_ticks()

    # model
    hess = np.log10(dataset['fit_hess'][plane_key])
    hess = np.ma.masked_invalid(hess, copy=True)
    im = ax_model.imshow(hess, vmin=0, vmax=3, **_imshow)
    ax_model.set_xlabel(x_label)
    for tl in ax_model.get_ymajorticklabels():
        tl.set_visible(False)
    cb = fig.colorbar(im, cax=ax_model_cb, orientation='horizontal')
    cb.set_label(r"$\log(N_*)$ Fit", bbox=label_bbox)
    cb.ax.xaxis.set_ticks_position('top')
    cb.locator = mpl.ticker.MultipleLocator(1.)
    cb.update_ticks()

    # chi
    hess = np.ma.masked_invalid(dataset['chi_hess'][plane_key], copy=True)
    im = ax_chi.imshow(hess, vmin=0, vmax=20, **_imshow)
    ax_chi.set_xlabel(x_label)
    for tl in ax_chi.get_ymajorticklabels():
        tl.set_visible(False)
    cb = fig.colorbar(im, cax=ax_chi_cb, orientation='horizontal')
    cb.set_label(r"$\chi^2$", bbox=label_bbox)
    cb.ax.xaxis.set_ticks_position('top')
    cb.locator = mpl.ticker.MultipleLocator(5)
    cb.update_ticks()

    # diff
    hess = np.ma.masked_invalid(dataset['diff_hess'][plane_key], copy=True)
    im = ax_diff.imshow(hess, vmin=-50, vmax=50, **_imshow_diff)
    for tl in ax_diff.get_ymajorticklabels():
        tl.set_visible(False)
    ax_diff.set_xlabel(x_label)
    cb = fig.colorbar(im, cax=ax_diff_cb, orientation='horizontal')
    cb.set_label(r"$\Delta_\mathrm{mock-model}$ ($N_*$)", bbox=label_bbox)
    cb.ax.xaxis.set_ticks_position('top')
    cb.locator = mpl.ticker.MultipleLocator(20)
    cb.update_ticks()

    # Tick hack. Ideally I'd want this data to be embedded in the the HDF5
    if plane_key == 'lewis':
        d_xticks = 0.5
        d_yticks = 1.
    else:  # OIR-ALL
        d_xticks = 2.
        d_yticks = 1.
    for i, ax in enumerate((ax_obs, ax_model, ax_chi, ax_diff)):
        ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(base=d_xticks))
        ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(base=d_yticks))
        ax.set_xlim(*x_span)
        ax.set_ylim(*y_span[::-1])
        if i > 0 and plane_key == 'lewis':
            ax.set_xticks([0.0, 0.5, 1.0])

    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def plot_sfh(model_sfh, mock_sfh, plot_path):
    labels = {'lewis': r'ACS-MS', 'oir_all': r'OIR-ALL'}
    colors = {'lewis': 'dodgerblue', 'oir_all': 'maroon'}

    fig = Figure(figsize=(3.5, 3.5), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(1, 1,
                           left=0.18, right=0.95, bottom=0.15, top=0.95,
                           wspace=None, hspace=None,
                           width_ratios=None, height_ratios=None)
    ax = fig.add_subplot(gs[0])
    for plane_key in ['lewis', 'oir_all']:
        if plane_key not in model_sfh['sfh'].keys():
            continue
        plot_single_sfh_line(ax, model_sfh['sfh'][plane_key],
                             label=labels[plane_key],
                             color=colors[plane_key],
                             drawstyle='steps-mid')
        _plot_mean_age(ax, model_sfh['sfh'][plane_key].attrs['mean_age'],
                       c=colors[plane_key])

    # plot_single_sfh_line(ax, model_sfh, label='Model', color='k')
    # print model_sfh['sfr']
    _plot_mock_sfh(ax, mock_sfh, lw=1.5, c='k', label='Mock')
    _plot_mean_age(ax, mock_sfh.attrs['mean_age'])

    ax.legend(loc='lower left', fontsize=8, frameon=True)

    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def _plot_mock_sfh(ax, table, **kwargs):
    A = table['log(age)']
    sfr = np.log10(table['sfr_msolar_yr'])
    _plot_args = {'drawstyle': 'steps-mid', 'label': 'Model'}
    _plot_args.update(**kwargs)
    ax.plot(A, sfr, **_plot_args)


def _plot_mean_age(ax, mean_age, **axvline_args):
    _axvline_args = {'lw': 0.8, 'ls': '--', 'c': 'k', 'zorder': 100}
    _axvline_args.update(axvline_args)
    ax.axvline(np.log10(mean_age[0]) + 9., **_axvline_args)


if __name__ == '__main__':
    main()
