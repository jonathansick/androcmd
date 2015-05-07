#!/usr/bin/env python
# encoding: utf-8
"""
Plotting utilities for androcmd
"""

import string
import numpy as np
from palettable.cubehelix import perceptual_rainbow_16
from palettable.colorbrewer.diverging import RdBu_11
import matplotlib as mpl
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.gridspec as gridspec
from astroML.plotting import scatter_contour


def contour_hess(ax, c, m, xlim, ylim,
                 threshold=20, levels=10, bins=100, log_counts=True,
                 plot_args=None, contour_args=None):
    """Plot a CMD as a contour Hess diagram in high-density regions, and
    as a scatter plot in low density regions.

    Parameters
    ----------
    ax :
        The matplotlib Axes instance.
    c : ndarray
        The colour (x) coordinates of stars
    m : ndarray
        The magnitude (y) coordinates of stars
    """
    default_plot_args = {'ms': 2.0, 'mfc': 'k', 'mec': 'None',
                         'rasterized': True, 'alpha': 0.3}
    if plot_args is not None:
        default_plot_args.update(plot_args)

    default_contour_args = {'cmap': mpl.cm.gray_r,
                            'linestyles': 'None',
                            'linewidths': 0.,
                            'alpha': 1.}
    if contour_args is not None:
        default_contour_args.append(contour_args)

    scatter_contour(c, m, levels=levels, threshold=threshold,
                    log_counts=log_counts,
                    histogram2d_args={'bins': bins,
                                      'range': [[min(xlim), max(xlim)],
                                                [min(ylim), max(ylim)]]},
                    plot_args=default_plot_args,
                    contour_args=default_contour_args,
                    ax=ax)


def plot_fit_grid(pipeline, dataset, fit_keys, plane_keys, plot_path,
                  ysize=3.5):
    n_y = len(fit_keys) + 1
    height_ratios = [0.1] + [1] * len(fit_keys)

    if len(fit_keys) > 1:
        multi_panel = True
    else:
        multi_panel = False

    fig = Figure(figsize=(7, ysize), frameon=False, dpi=300)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(n_y, 4, wspace=0.15, hspace=0.2,
                           left=0.08, bottom=0.15, right=0.95,
                           width_ratios=(1, 1, 1, 1),
                           height_ratios=height_ratios)

    for i, (fit_key, plane_key) in enumerate(zip(fit_keys, plane_keys)):
        if i == n_y - 2:
            last = True
        else:
            last = False
        _plot_plane(pipeline, dataset, fit_key, plane_key, i, fig, gs,
                    last=last, multi_panel=multi_panel)

    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def _plot_plane(pipeline, dataset, fit_key, plane_key, i, fig, gs,
                last=False, multi_panel=True):
    obs_hess = pipeline.make_obs_hess(dataset, plane_key)
    fit_hess = pipeline.make_fit_hess(fit_key, plane_key)
    sigma = np.sqrt(obs_hess.hess)
    chi = ((obs_hess.hess - fit_hess.hess) / sigma) ** 2.
    diff = obs_hess.hess - fit_hess.hess

    ax_obs = fig.add_subplot(gs[i + 1, 0])
    ax_model = fig.add_subplot(gs[i + 1, 1])
    ax_chi = fig.add_subplot(gs[i + 1, 2])
    ax_diff = fig.add_subplot(gs[i + 1, 3])

    cube_map = perceptual_rainbow_16.mpl_colormap
    div_map = RdBu_11.mpl_colormap

    fit_map = pipeline.plot_fit_hess(ax_model, fit_key, plane_key,
                                     imshow=dict(vmin=0, vmax=3.,
                                                 cmap=cube_map))
    ax_model.yaxis.set_major_formatter(mpl.ticker.NullFormatter())
    ax_model.set_ylabel('')

    obs_map = pipeline.plot_obs_hess(ax_obs, dataset, plane_key,
                                     imshow=dict(vmin=0, vmax=3.,
                                                 cmap=cube_map))

    chi_map = pipeline.plot_hess_array(ax_chi, chi, plane_key, log=False,
                                       imshow=dict(vmax=20, cmap=cube_map))
    ax_chi.yaxis.set_major_formatter(mpl.ticker.NullFormatter())
    ax_chi.set_ylabel('')

    diff_map = pipeline.plot_hess_array(ax_diff, diff, plane_key, log=False,
                                        imshow=dict(vmin=-50, vmax=50,
                                                    cmap=div_map))
    ax_diff.yaxis.set_major_formatter(mpl.ticker.NullFormatter())
    ax_diff.set_ylabel('')

    if not last:
        ax_diff.set_xlabel('')
        ax_chi.set_xlabel('')
        ax_model.set_xlabel('')
        ax_obs.set_xlabel('')
        ax_obs.xaxis.set_major_formatter(mpl.ticker.NullFormatter())
        ax_model.xaxis.set_major_formatter(mpl.ticker.NullFormatter())
        ax_chi.xaxis.set_major_formatter(mpl.ticker.NullFormatter())
        ax_diff.xaxis.set_major_formatter(mpl.ticker.NullFormatter())

    if i == 0:  # colorbar for first row only
        ax_obs_cb = fig.add_subplot(gs[0, 0])
        ax_model_cb = fig.add_subplot(gs[0, 1])
        ax_chi_cb = fig.add_subplot(gs[0, 2])
        ax_diff_cb = fig.add_subplot(gs[0, 3])

        obs_cb = fig.colorbar(obs_map, cax=ax_obs_cb, orientation='horizontal')
        obs_cb.set_label(r"$\log(N_*)$ Obs.", size=9)
        obs_cb.ax.xaxis.set_ticks_position('top')
        obs_cb.locator = mpl.ticker.MultipleLocator(1.0)
        for tl in obs_cb.ax.get_xmajorticklabels():
            tl.set_size(8.)
        obs_cb.update_ticks()

        fit_cb = fig.colorbar(fit_map, cax=ax_model_cb,
                              orientation='horizontal')
        fit_cb.set_label(r"$\log(N_*)$ Model", size=9)
        fit_cb.ax.xaxis.set_ticks_position('top')
        fit_cb.locator = mpl.ticker.MultipleLocator(1.0)
        for tl in fit_cb.ax.get_xmajorticklabels():
            tl.set_size(8.)
        fit_cb.update_ticks()

        chi_cb = fig.colorbar(chi_map, cax=ax_chi_cb, orientation='horizontal')
        chi_cb.set_label(r"$\chi^2$", size=9)
        chi_cb.ax.xaxis.set_ticks_position('top')
        chi_cb.locator = mpl.ticker.MultipleLocator(5)
        for tl in chi_cb.ax.get_xmajorticklabels():
            tl.set_size(8.)
        chi_cb.update_ticks()

        diff_cb = fig.colorbar(diff_map, cax=ax_diff_cb,
                               orientation='horizontal')
        diff_cb.set_label(r"$\Delta_\mathrm{obs-model}$ ($N_*$)", size=9)
        diff_cb.ax.xaxis.set_ticks_position('top')
        diff_cb.locator = mpl.ticker.MultipleLocator(20)
        for tl in diff_cb.ax.get_xmajorticklabels():
            tl.set_size(8.)
        diff_cb.update_ticks()

    if multi_panel:
        # more than one row; add subfig annotations
        alphanum = dict(zip(range(1, 27), string.ascii_lowercase))
        alpha = alphanum[i + 1]
        txt = '({0})'.format(alpha)
        ax_obs.text(-0.38, 1.0, txt,
                    transform=ax_obs.transAxes,
                    ha='left',
                    va='top',
                    size=11)
