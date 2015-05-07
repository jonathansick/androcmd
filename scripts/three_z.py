#!/usr/bin/env python
# encoding: utf-8
"""
Modelling with three metallicity groups.
"""

import argparse

from starfisher.pipeline import PipelineBase
from androcmd.planes import BasicPhatPlanes
from androcmd.phatpipeline import (
    ExtendedSolarIsocs, ExtendedSolarLockfile,
    NoDust, PhatCrowding)
from androcmd.phatpipeline import PhatCatalog


def main():
    args = parse_args()

    pipeline = Pipeline(root_dir="b{0:d}_threez".format(args.brick),
                        isoc_args=dict(isoc_kind='parsec_CAF09_v1.2S',
                                       photsys_version='yang'))
    fit_planes = {'f475w_f160w': ['f475w_f160w'],
                  'rgb': ['f475w_f814w_rgb'],
                  'ms': ['f475w_f814w_ms']}
    if len(fit_planes) > 0:
        dataset = PhatCatalog(args.brick)
    for fit_name in args.fit:
        plane = fit_planes[fit_name]
        pipeline.fit(fit_name, plane, dataset)

    if args.plot:
        plot_fit(pipeline, dataset)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Multi-metallicity fitting")
    parser.add_argument('brick', type=int)
    parser.add_argument('--fit', nargs='*', default=list())
    parser.add_argument('--plot', default=False, action='store_true')
    return parser.parse_args()


class Pipeline(BasicPhatPlanes,
               ExtendedSolarIsocs, ExtendedSolarLockfile,
               NoDust, PhatCrowding, PipelineBase):
    """A pipeline for fitting PHAT bricks with solar metallicity isochrones."""
    def __init__(self, **kwargs):
        print "MultiZPhatPipeline", kwargs
        super(Pipeline, self).__init__(**kwargs)


def plot_fit(pipeline, dataset):
    """Plot for the thesis"""
    fit_key = 'f475w_f160w'
    plane_key = 'f475w_f160w'
    plotpath = "b23_three_z_f475w_f160w"

    import numpy as np
    from palettable.cubehelix import perceptual_rainbow_16
    from palettable.colorbrewer.diverging import RdBu_11
    import matplotlib as mpl
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    import matplotlib.gridspec as gridspec

    cube_map = perceptual_rainbow_16.mpl_colormap

    obs_hess = pipeline.make_obs_hess(dataset, plane_key)
    fit_hess = pipeline.make_fit_hess(fit_key, plane_key)
    sigma = np.sqrt(obs_hess.hess)
    chi = ((obs_hess.hess - fit_hess.hess) / sigma) ** 2.
    diff = obs_hess.hess - fit_hess.hess

    fig = Figure(figsize=(7, 3.5), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(2, 4, wspace=0.15, hspace=0.2,
                           left=0.08, bottom=0.15, right=0.95,
                           width_ratios=(1, 1, 1, 1), height_ratios=(0.1, 1))
    ax_obs = fig.add_subplot(gs[1, 0])
    ax_model = fig.add_subplot(gs[1, 1])
    ax_chi = fig.add_subplot(gs[1, 2])
    ax_diff = fig.add_subplot(gs[1, 3])
    ax_obs_cb = fig.add_subplot(gs[0, 0])
    ax_model_cb = fig.add_subplot(gs[0, 1])
    ax_chi_cb = fig.add_subplot(gs[0, 2])
    ax_diff_cb = fig.add_subplot(gs[0, 3])

    fit_map = pipeline.plot_fit_hess(ax_model, fit_key, plane_key,
                                     imshow=dict(vmin=0, vmax=3.,
                                                 cmap=cube_map))
    ax_model.yaxis.set_major_formatter(mpl.ticker.NullFormatter())
    ax_model.set_ylabel('')
    fit_cb = fig.colorbar(fit_map, cax=ax_model_cb,
                          orientation='horizontal')
    fit_cb.set_label(r"$\log(N_*)$ Model", size=9)
    fit_cb.ax.xaxis.set_ticks_position('top')
    fit_cb.locator = mpl.ticker.MultipleLocator(1.0)
    for tl in fit_cb.ax.get_xmajorticklabels():
        tl.set_size(8.)
    fit_cb.update_ticks()

    obs_map = pipeline.plot_obs_hess(ax_obs, dataset, plane_key,
                                     imshow=dict(vmin=0, vmax=3.,
                                                 cmap=cube_map))
    ax_obs.yaxis.set_major_formatter(mpl.ticker.NullFormatter())
    obs_cb = fig.colorbar(obs_map, cax=ax_obs_cb, orientation='horizontal')
    obs_cb.set_label(r"$\log(N_*)$ Obs.", size=9)
    obs_cb.ax.xaxis.set_ticks_position('top')
    obs_cb.locator = mpl.ticker.MultipleLocator(1.0)
    for tl in obs_cb.ax.get_xmajorticklabels():
        tl.set_size(8.)
    obs_cb.update_ticks()

    chi_map = pipeline.plot_hess_array(ax_chi, chi, plane_key, log=False,
                                       imshow=dict(vmax=20, cmap=cube_map))
    ax_chi.yaxis.set_major_formatter(mpl.ticker.NullFormatter())
    ax_chi.set_ylabel('')
    chi_cb = fig.colorbar(chi_map, cax=ax_chi_cb, orientation='horizontal')
    chi_cb.set_label(r"$\chi^2$", size=9)
    chi_cb.ax.xaxis.set_ticks_position('top')
    chi_cb.locator = mpl.ticker.MultipleLocator(5)
    for tl in chi_cb.ax.get_xmajorticklabels():
        tl.set_size(8.)
    chi_cb.update_ticks()

    diff_map = pipeline.plot_hess_array(ax_diff, diff, plane_key, log=False,
                                        imshow=dict(vmin=-50, vmax=50,
                                                    cmap=RdBu_11.mpl_colormap))
    ax_diff.yaxis.set_major_formatter(mpl.ticker.NullFormatter())
    ax_diff.set_ylabel('')
    diff_cb = fig.colorbar(diff_map, cax=ax_diff_cb,
                           orientation='horizontal')
    diff_cb.set_label(r"$\Delta_\mathrm{obs-model}$ ($N_*$)", size=9)
    diff_cb.ax.xaxis.set_ticks_position('top')
    diff_cb.locator = mpl.ticker.MultipleLocator(20)
    for tl in diff_cb.ax.get_xmajorticklabels():
        tl.set_size(8.)
    diff_cb.update_ticks()

    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plotpath + ".pdf", format="pdf")


if __name__ == '__main__':
    main()
