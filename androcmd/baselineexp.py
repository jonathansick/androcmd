# encoding: utf-8
"""
Helper code for the baseline experiment.

2015-05-15 - Created by Jonathan Sick
"""

import os
from collections import OrderedDict

import numpy as np

import matplotlib as mpl
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.gridspec as gridspec
from palettable.cubehelix import perceptual_rainbow_16
from palettable.colorbrewer.diverging import RdBu_11
import palettable

from starfisher.pipeline import PipelineBase
from starfisher.sfhplot import plot_single_sfh_line

from androcmd.planes import BaselineTestPhatPlanes
from androcmd.phatpipeline import (SolarZIsocs, SolarLockfile,
                                   LewisBrickDust, PhatCrowding,
                                   ExtendedSolarIsocs, ExtendedSolarLockfile)


class SolarZPipeline(BaselineTestPhatPlanes, SolarZIsocs,
                     SolarLockfile, LewisBrickDust, PhatCrowding,
                     PipelineBase):
    """Pipeline for the baseline test with only solar metallicity."""
    def __init__(self, **kwargs):
        super(SolarZPipeline, self).__init__(**kwargs)


class ThreeZPipeline(BaselineTestPhatPlanes, ExtendedSolarIsocs,
                     ExtendedSolarLockfile, LewisBrickDust, PhatCrowding,
                     PipelineBase):
    """Pipeline for baseline test with three metallicity tracks."""
    def __init__(self, **kwargs):
        super(ThreeZPipeline, self).__init__(**kwargs)


def init_fits(p, fit_labels, dataset):
    base_dir = os.path.join(os.getenv('STARFISH'), p.root_dir)
    useable_fits = []
    for label in fit_labels:
        fit_dir = os.path.join(base_dir, label)
        print "fit_dir", fit_dir
        if os.path.exists(fit_dir):
            useable_fits.append(label)
            if label == 'lewis':
                p.fit('lewis', ['lewis'], dataset)
            elif label == 'acs_rgb':
                p.fit('acs_rgb', ['acs_rgb'], dataset)
            elif label == 'acs_all':
                p.fit('acs_all', ['acs_all'], dataset)
            elif label == 'oir_all':
                p.fit('oir_all', ['oir_all'], dataset)
            elif label == 'ir_rgb':
                p.fit('ir_rgb', ['ir_rgb'], dataset)
    return useable_fits


def init_grid_plot(plot_path, p, dataset):
    fit_labels = OrderedDict((
        ('lewis', 'Fitting ACS-MS'),
        ('acs_rgb', 'Fitting ACS-RGB'),
        ('acs_all', 'Fitting ACS-ALL'),
        ('oir_all', 'Fitting OIR-ALL'),
        ('ir_rgb', 'Fitting NIR-RGB')))
    nfits = len(fit_labels)
    nplanes = len(fit_labels)

    xlocators = {'lewis': mpl.ticker.MultipleLocator(base=0.4),
                 'acs_rgb': mpl.ticker.MultipleLocator(base=1.),
                 'acs_all': mpl.ticker.MultipleLocator(base=2.),
                 'oir_all': mpl.ticker.MultipleLocator(base=2.),
                 'ir_rgb': mpl.ticker.MultipleLocator(base=0.3)}

    # only re-initialize fits that are available
    useable_fits = init_fits(p, fit_labels, dataset)

    fig = Figure(figsize=(7.5, 8.5), frameon=False)
    canvas = FigureCanvas(fig)
    # TODO add a colorbar axis
    gs = gridspec.GridSpec(nfits, nplanes + 1,
                           left=0.04, right=0.92, bottom=0.05, top=0.97,
                           wspace=0.3, hspace=0.25,
                           width_ratios=[1.] * nplanes + [0.1],
                           height_ratios=None)

    cb_ax = fig.add_subplot(gs[:, -1])
    axes = {}
    for i, fit_key in enumerate(fit_labels):
        fit_axes = {k: fig.add_subplot(gs[i, j])
                    for j, k in enumerate(fit_labels)}
        axes[fit_key] = fit_axes

    # Ensure even blank axes have consistent extents+labels
    for fit_key in fit_labels:
        for plane_key in fit_labels:
            ax = axes[fit_key][plane_key]
            p.init_plane_axes(ax, plane_key)

    # x-labels only for the bottom row
    for fit_key in fit_labels.keys()[:-1]:
        for plane_key in fit_labels.keys():
            ax = axes[fit_key][plane_key]
            for tl in ax.get_xmajorticklabels():
                tl.set_visible(False)
            ax.set_xlabel('')

    # Make smaller tick labels
    for fit_key in fit_labels.keys():
        for plane_key in fit_labels.keys():
            ax.tick_params(axis='both', which='major', labelsize=7)
            ax.tick_params(axis='both', which='minor', labelsize=7)
            # Move the y-label into the plot
            ax.yaxis.labelpad = -100
            ax.xaxis.label.set_size(7)
            ax.yaxis.label.set_size(7)

    # Label the fit name in each row
    for fit_key in fit_labels.keys():
        plane_key = fit_labels.keys()[0]
        ax = axes[fit_key][plane_key]
        ax.text(0.00, 1.05, fit_labels[fit_key],
                ha='left', va='baseline',
                transform=ax.transAxes,
                size=10)

    # Highlight the fitted Hess plane in each row
    for fit_key in fit_labels:
        for plane_key in fit_labels:
            ax = axes[fit_key][plane_key]

            if fit_key == plane_key:
                highlight_color = '#3498db'
                for loc in ['bottom', 'top', 'right', 'left']:
                    ax.spines[loc].set_color(highlight_color)
                    ax.spines[loc].set_linewidth(3.5)

            ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(base=1))
            ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(base=0.25))
            ax.xaxis.set_major_locator(xlocators[plane_key])
            ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(base=0.1))

    return fit_labels, useable_fits, \
        fig, canvas, axes, cb_ax


def reset_grid_plot_labels(fit_labels, axes):
    # x-labels only for the bottom row
    for fit_key in fit_labels.keys()[:-1]:
        for plane_key in fit_labels.keys():
            ax = axes[fit_key][plane_key]
            for tl in ax.get_xmajorticklabels():
                tl.set_visible(False)
            ax.set_xlabel('')

    # Make smaller tick labels
    for fit_key in fit_labels.keys():
        for plane_key in fit_labels.keys():
            ax = axes[fit_key][plane_key]
            ax.tick_params(axis='both', which='major', labelsize=7)
            ax.tick_params(axis='both', which='minor', labelsize=7)
            # Move the y-label into the plot
            ax.yaxis.labelpad = -26
            ax.xaxis.label.set_size(9)
            ax.yaxis.label.set_size(9)


def plot_fit_hess_grid(plot_path, p, dataset):
    fit_labels, useable_fits, \
        fig, canvas, axes, cb_ax = init_grid_plot(plot_path, p, dataset)

    cube_map = perceptual_rainbow_16.mpl_colormap
    imshow_args = dict(vmax=20, cmap=cube_map)

    for fit_key in useable_fits:
        for plane_key in fit_labels:
            chi_hess = p.make_chisq_hess(dataset, fit_key, plane_key)
            chi_data = chi_hess.masked_hess
            chi_red = p.compute_fit_chi(dataset, fit_key, plane_key,
                                        chi_hess=chi_hess)
            print fit_key, plane_key, chi_red

            ax = axes[fit_key][plane_key]
            chi_map = p.plot_hess_array(ax,
                                        chi_data,
                                        plane_key,
                                        log=False,
                                        imshow=imshow_args)

            if chi_red < 100.:
                txt = '$\chi^2_\mathrm{{r}}={0:.1f}$'.format(chi_red)
            else:
                txt = '$\chi^2_\mathrm{{r}}={0:.1e}$'.format(chi_red)
                a, b = txt.split('e+')
                b = b.rstrip('$')
                b = int(b)
                txt = '{0} \\times 10^{{{1:d}}}$'.format(a, b)
            ax.text(0.05, 0.95,
                    txt,
                    ha='left', va='top',
                    backgroundcolor='w',
                    size=9,
                    transform=ax.transAxes)

    cb = fig.colorbar(ax=ax, cax=cb_ax, mappable=chi_map)
    cb.set_label(r'$\chi^2$')

    reset_grid_plot_labels(fit_labels, axes)

    canvas.print_figure(plot_path + ".pdf", format="pdf")


def plot_diff_hess_grid(plot_path, p, dataset):
    fit_labels, useable_fits, \
        fig, canvas, axes, cb_ax = init_grid_plot(plot_path, p, dataset)

    div_map = RdBu_11.mpl_colormap
    imshow_args = dict(vmin=-50, vmax=50, cmap=div_map)

    for fit_key in useable_fits:
        for plane_key in fit_labels:
            chi_hess = p.make_chisq_hess(dataset, fit_key, plane_key)
            chi_red = p.compute_fit_chi(dataset, fit_key, plane_key,
                                        chi_hess=chi_hess)
            print fit_key, plane_key, chi_red

            # compute a difference map
            obs_hess = p.make_obs_hess(dataset, plane_key)
            fit_hess = p.make_fit_hess(fit_key, plane_key)
            diff = obs_hess.masked_hess - fit_hess.masked_hess

            ax = axes[fit_key][plane_key]
            chi_map = p.plot_hess_array(ax,
                                        diff,
                                        plane_key,
                                        log=False,
                                        imshow=imshow_args)

            if chi_red < 100.:
                txt = '$\chi^2_\mathrm{{r}}={0:.1f}$'.format(chi_red)
            else:
                txt = '$\chi^2_\mathrm{{r}}={0:.1e}$'.format(chi_red)
                a, b = txt.split('e+')
                b = b.rstrip('$')
                b = int(b)
                txt = '{0} \\times 10^{{{1:d}}}$'.format(a, b)
            ax.text(0.05, 0.95,
                    txt,
                    ha='left', va='top',
                    backgroundcolor='w',
                    size=9,
                    transform=ax.transAxes)

    cb = fig.colorbar(ax=ax, cax=cb_ax, mappable=chi_map)
    cb.set_label(r"$\Delta_\mathrm{obs-model}$ ($N_*$)")

    reset_grid_plot_labels(fit_labels, axes)

    canvas.print_figure(plot_path + ".pdf", format="pdf")


def tabulate_fit_chi(table_path, p, dataset):
    pass

    # fit_labels = OrderedDict((
    #     ('lewis', 'Fitting ACS-MS'),
    #     ('acs_rgb', 'Fitting ACS-RGB'),
    #     ('acs_all', 'Fitting ACS-ALL'),
    #     ('oir_all', 'Fitting OIR-ALL'),
    #     ('ir_rgb', 'Fitting NIR-RGB')))

    # only re-initialize fits that are available
    # useable_fits = init_fits(p, fit_labels)

    # nfits = len(fit_labels)
    # nplanes = len(fit_labels)


def sfh_comparison_plot(plot_path, p, dataset):
    fit_labels = OrderedDict((
        ('lewis', 'Fitting ACS-MS'),
        ('acs_rgb', 'Fitting ACS-RGB'),
        ('acs_all', 'Fitting ACS-ALL'),
        ('oir_all', 'Fitting OIR-ALL'),
        ('ir_rgb', 'Fitting NIR-RGB')))

    # only re-initialize fits that are available
    useable_fits = init_fits(p, fit_labels, dataset)

    fig = Figure(figsize=(7.5, 3.5), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(1, 1,
                           left=0.08, right=0.95, bottom=0.15, top=0.95,
                           wspace=None, hspace=None,
                           width_ratios=None, height_ratios=None)
    ax = fig.add_subplot(gs[0])
    colors = dict(zip(fit_labels.keys(),
                      palettable.tableau.ColorBlind_10.mpl_colors))
    hatches = dict(zip(fit_labels.keys(),
                       ['//', '\\\\', '||', '--', 'oo']))
    for fit_key in useable_fits:
        print "fit_key", fit_key
        sfh_table = p.fits[fit_key].solution_table(marginalize_z=True)
        mean_age, mean_age_sigma = p.fits[fit_key].mean_age
        # mean_log_age = p.fits[fit_key].mean_log_age
        plot_single_sfh_line(
            ax, sfh_table,
            z_formatter=mpl.ticker.FormatStrFormatter("%.2f"),
            age_formatter=mpl.ticker.FormatStrFormatter("%4.1f"),
            color=colors[fit_key],
            label=r'{0} $\langle A \rangle={1:.1f}\pm{2:.1f}$ Gyr'.format(
                fit_labels[fit_key], mean_age, mean_age_sigma),
            age_lim=(1e-3, 14.),
            amp_key='sfr',
            log_amp=True,
            log_age=True,
            x_label=True,
            y_label=True,
            plot_errors=True,
            hatch_errors=hatches[fit_key])

    ax.legend(frameon=True, loc='lower left', fontsize=8,
              fancybox=True, framealpha=0.5)

    for logage in np.log10(np.arange(1, 14, 1) * 1e9):
        ax.axvline(logage, c='r', ls='-', lw=0.7, zorder=-20)

    ax.set_ylim(-9, 6.)
    ax.set_xlim(6.5, 10.2)

    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")
