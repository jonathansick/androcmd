#!/usr/bin/env python
# encoding: utf-8
'''
Plot mock tests made by the mock_sfh_test.py script.

2015-07-09 - Created by Jonathan Sick
'''

import os
import argparse

import numpy as np
import h5py

import matplotlib as mpl
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.gridspec as gridspec

from palettable.cubehelix import perceptual_rainbow_16

from starfisher.sfhplot import plot_single_sfh_line


def main():
    args = parse_args()
    dirname = os.path.join(os.getenv('STARFISH'), args.name)
    h5path = os.path.join(dirname, args.name + '.hdf5')
    hdf5 = h5py.File(h5path, 'r')

    plot_hess_planes(hdf5['mocksfh'], dirname)
    plot_sfhs(hdf5['mocksfh'], dirname)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('name',
                        help='Test name')
    return parser.parse_args()


def plot_sfhs(dataset, base_dir):
    sfh_list = dataset.keys()
    for sfh_run in sfh_list:
        plot_path = os.path.join(base_dir,
                                 '{0}_sfh'.format(sfh_run))
        plot_sfh(dataset[sfh_run], dataset[sfh_run]['model_sfh_marginal'],
                 plot_path)


def plot_hess_planes(dataset, base_dir):
    sfh_list = dataset.keys()
    for sfh_run in sfh_list:
        plane_keys = dataset[sfh_run]['sim_hess'].keys()
        for plane_key in plane_keys:
            plot_path = os.path.join(base_dir,
                                     '{0}_{1}_hess'.format(sfh_run, plane_key))
            _plot_hess(dataset[sfh_run], plane_key, plot_path)


def _plot_hess(dataset, plane_key, plot_path):
    fig = Figure(figsize=(6.5, 3.5), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(2, 4,
                           left=0.07, right=0.95, bottom=0.15, top=0.90,
                           wspace=None, hspace=0.,
                           width_ratios=None, height_ratios=(0.1, 1))

    extent = dataset['sim_hess'][plane_key].attrs['extent']
    origin = dataset['sim_hess'][plane_key].attrs['origin']
    y_label = dataset['sim_hess'][plane_key].attrs['y_label']
    x_label = dataset['sim_hess'][plane_key].attrs['x_label']

    cube_map = perceptual_rainbow_16.mpl_colormap

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

    # observations
    hess = np.log10(dataset['obs_hess'][plane_key])
    hess = np.ma.masked_invalid(hess, copy=True)
    im = ax_obs.imshow(hess, vmin=0, vmax=3, **_imshow)
    ax_obs.set_xlabel(x_label)
    ax_obs.set_ylabel(y_label)
    cb = fig.colorbar(im, cax=ax_obs_cb, orientation='horizontal')
    cb.set_label(r"$\log(N_*)$ Obs.")
    cb.ax.xaxis.set_ticks_position('top')
    cb.locator = mpl.ticker.MultipleLocator(1.)
    cb.update_ticks()

    # model
    h = np.array(dataset['sim_hess'][plane_key][:, :])
    print "sim_hess", plot_path, np.nanmin(h), np.nanmax(h)
    hess = np.log10(dataset['sim_hess'][plane_key])
    hess = np.ma.masked_invalid(hess, copy=True)
    im = ax_model.imshow(hess, vmin=-6, vmax=-3, **_imshow)
    ax_model.set_xlabel(x_label)
    for tl in ax_model.get_ymajorticklabels():
        tl.set_visible(False)
    cb = fig.colorbar(im, cax=ax_model_cb, orientation='horizontal')
    cb.set_label(r"$\log(N_*)$ Model")
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
    cb.set_label(r"$\chi^2$")
    cb.ax.xaxis.set_ticks_position('top')
    cb.locator = mpl.ticker.MultipleLocator(5)
    cb.update_ticks()

    # diff
    hess = np.ma.masked_invalid(dataset['diff_hess'][plane_key], copy=True)
    im = ax_diff.imshow(hess, vmin=-50, vmax=50, **_imshow)
    for tl in ax_diff.get_ymajorticklabels():
        tl.set_visible(False)
    ax_diff.set_xlabel(x_label)
    cb = fig.colorbar(im, cax=ax_diff_cb, orientation='horizontal')
    cb.set_label(r"$\Delta_\mathrm{obs-model}$ ($N_*$)")
    cb.ax.xaxis.set_ticks_position('top')
    cb.locator = mpl.ticker.MultipleLocator(20)
    cb.update_ticks()

    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def plot_sfh(mock_sfh, model_sfh, plot_path):
    labels = {'lewis': r'ACS-MS', 'oir_all': r'OIR-ALL'}
    colors = {'lewis': 'dodgerblue', 'oir_all': 'maroon'}

    fig = Figure(figsize=(3.5, 3.5), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(1, 1,
                           left=0.15, right=0.95, bottom=0.15, top=0.95,
                           wspace=None, hspace=None,
                           width_ratios=None, height_ratios=None)
    ax = fig.add_subplot(gs[0])
    for plane_key in ['lewis', 'oir_all']:
        plot_single_sfh_line(ax, mock_sfh['sfh'][plane_key],
                             label=labels[plane_key],
                             color=colors[plane_key])

    # plot_single_sfh_line(ax, model_sfh, label='Model', color='k')
    # print model_sfh['sfr']
    _plot_mock_sfh(ax, model_sfh, lw=1.5, c='k', label='Model')

    ax.legend(loc='lower right', fontsize=8, frameon=True)

    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def _plot_mock_sfh(ax, table, **kwargs):
    A = table['log(age)']
    sfr = np.log10(table['sfr_msolar_yr'])
    _plot_args = {'drawstyle': 'steps-mid', 'label': 'Model'}
    _plot_args.update(**kwargs)
    ax.plot(A, sfr, **_plot_args)


if __name__ == '__main__':
    main()
