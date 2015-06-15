# encoding: utf-8
"""
Helper code for the baseline experiment.

2015-05-15 - Created by Jonathan Sick
"""

import os
from collections import OrderedDict

import numpy as np
from astropy.coordinates import Distance
import astropy.units as u

import matplotlib as mpl
import matplotlib.pyplot as plt
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
from androcmd.phatpipeline import get_demo_age_grid


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
                      palettable.colorbrewer.qualitative.Set1_5.mpl_colors))
    # hatches = dict(zip(fit_labels.keys(),
    #                    ['//', '\\\\', '||', '--', 'oo']))
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
            hatch_errors=None)

    ax.legend(frameon=True, loc='lower center', fontsize=8,
              fancybox=True, framealpha=0.5)

    for logage in np.log10(np.arange(1, 14, 1) * 1e9):
        ax.axvline(logage, c='r', ls='-', lw=0.7, zorder=-20)

    ax.set_ylim(-9, 5.)
    ax.set_xlim(6.5, 10.2)

    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def plot_sfh_metallicity_trends(plot_path, p, dataset, fit_key):
    fit_labels = OrderedDict((
        ('lewis', 'Fitting ACS-MS'),
        ('acs_rgb', 'Fitting ACS-RGB'),
        ('acs_all', 'Fitting ACS-ALL'),
        ('oir_all', 'Fitting OIR-ALL'),
        ('ir_rgb', 'Fitting NIR-RGB')))

    # only re-initialize fits that are available
    init_fits(p, fit_labels, dataset)

    fig = Figure(figsize=(3.5, 3.5), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(1, 1,
                           left=0.15, right=0.95, bottom=0.15, top=0.95,
                           wspace=None, hspace=None,
                           width_ratios=None, height_ratios=None)
    ax = fig.add_subplot(gs[0])

    colors = palettable.tableau.ColorBlind_10.mpl_colors

    # TODO implement these
    sfh_tables = p.fits[fit_key].solution_table(split_z=True)
    mean_ages, mean_age_sigmas = p.fits[fit_key].mean_age_by_z
    legend_str = r'$\log Z/Z_\odot = {0}$, ' \
        r'$\langle A \rangle={1:.1f}\pm{2:.1f}$ Gyr'
    for (z_key, color) in zip(sfh_tables.keys(), colors):
        sfh_table = sfh_tables[z_key]
        mean_age = mean_ages[z_key]
        mean_age_sigma = mean_age_sigmas[z_key]
        plot_single_sfh_line(
            ax, sfh_table,
            z_formatter=mpl.ticker.FormatStrFormatter("%.2f"),
            age_formatter=mpl.ticker.FormatStrFormatter("%4.1f"),
            color=color,
            label=legend_str.format(z_key, mean_age, mean_age_sigma),
            age_lim=(1e-3, 14.),
            amp_key='sfr',
            log_amp=True,
            log_age=True,
            x_label=True,
            y_label=True,
            plot_errors=True,
            hatch_errors=None)

    ax.legend(frameon=True, loc='best', fontsize=8,
              fancybox=True, framealpha=0.8)

    for logage in np.log10(np.arange(1, 14, 1) * 1e9):
        ax.axvline(logage, c='r', ls='-', lw=0.7, zorder=-20)

    ax.set_ylim(-9, 5.)
    ax.set_xlim(6.5, 10.2)

    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def plot_isocs(plot_path, pipeline, dataset):
    fig = Figure(figsize=(6.5, 5.), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(2, 3,
                           left=0.08, right=0.85, bottom=0.08, top=0.95,
                           wspace=0.15, hspace=0.25,
                           width_ratios=(1, 1, 0.1), height_ratios=(0.1, 1.))
    cax_ages = fig.add_subplot(gs[0, 0])
    cax_phases = fig.add_subplot(gs[1, 2])
    ax_ages = fig.add_subplot(gs[1, 0])
    ax_phases = fig.add_subplot(gs[1, 1])

    isoc_set = get_demo_age_grid(**dict(isoc_kind='parsec_CAF09_v1.2S',
                                        photsys_version='yang'))
    plane_key = 'oir_all'
    # plane = pipeline.planes[plane_key]

    # Plot the observed Hess diagram  in each axes
    for ax in [ax_ages, ax_phases]:
        pipeline.plot_obs_hess(ax_ages, dataset, plane_key, imshow=None)
        pipeline.plot_obs_hess(ax_phases, dataset, plane_key, imshow=None)

    # Plot isochrones by age
    cmap = palettable.cubehelix.perceptual_rainbow_16.mpl_colormap
    scalar_map = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=7.,
                                       vmax=10.1),
                                       cmap=cmap)
    scalar_map.set_array(np.array([isoc.age for isoc in isoc_set]))

    d = Distance(785 * u.kpc)
    for isoc in isoc_set:
        ax_ages.plot(isoc['F475W'] - isoc['F160W'],
                     isoc['F160W'] + d.distmod.value,
                     c=scalar_map.to_rgba(np.log10(isoc.age)))
    cax_ages = plt.colorbar(mappable=scalar_map, cax=cax_ages, ax=ax_ages,
                            orientation='horizontal')
    cax_ages.set_label(r"$\log(A/\mathrm{yr})$")

    # Plot phases
    phase_labels = {0: 'Pre-MS', 1: 'MS', 2: 'SGB', 3: 'RGB',
                    4: 'CHeB(1)', 5: 'CHeB(2)', 6: 'CHeB(3)',
                    7: 'E-AGB', 8: 'TP-AGB'}
    cmap = mpl.colors.ListedColormap(
        palettable.colorbrewer.qualitative.Set1_9.mpl_colors)
    scalar_map = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=-0.5,
                                                                 vmax=8.5),
                                       cmap=cmap)
    scalar_map.set_array(np.array(range(0, 9)))

    d = Distance(785 * u.kpc)
    for isoc in isoc_set:
        phases = np.unique(isoc['stage'])
        srt = np.argsort(phases)
        phases = phases[srt]
        for p in phases:
            s = np.where(isoc['stage'] == p)[0]
            ax_phases.plot(isoc['F475W'][s] - isoc['F160W'][s],
                           isoc['F160W'][s] + d.distmod.value,
                           c=scalar_map.to_rgba(p),
                           lw=0.8)
    cb_phases = plt.colorbar(mappable=scalar_map,
                             cax=cax_phases, ax=ax_phases, ticks=range(0, 9),
                             orientation='vertical')
    # for tl in ax.get_xmajorticklabels():
    #     tl.set_visible(False)
    # for label in cb_phases.ax.get_xmajorticklabels():
    #     label.set_rotation('vertical')
    cb_phases.ax.set_yticklabels([phase_labels[p] for p in range(0, 9)])
    cb_phases.set_label(r"Stage")
    # cb_phases.update_ticks()

    for tl in ax_phases.get_ymajorticklabels():
        tl.set_visible(False)
    ax_phases.set_ylabel('')

    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")
