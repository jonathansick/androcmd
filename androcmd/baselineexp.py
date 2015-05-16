# encoding: utf-8
"""
Helper code for the baseline experiment.

2015-05-15 - Created by Jonathan Sick
"""

from collections import OrderedDict

import matplotlib as mpl
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.gridspec as gridspec
from palettable.cubehelix import perceptual_rainbow_16

from starfisher.pipeline import PipelineBase
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


def plot_fit_hess_grid(plot_path, p, dataset):
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

    # p.fit('lewis', ['lewis'], dataset)
    p.fit('acs_rgb', ['acs_rgb'], dataset)
    p.fit('acs_all', ['acs_all'], dataset)
    p.fit('oir_all', ['oir_all'], dataset)
    p.fit('ir_rgb', ['ir_rgb'], dataset)

    usable_fits = ['acs_rgb', 'acs_all', 'ir_rgb', 'oir_all']

    fig = Figure(figsize=(7.5, 8.5), frameon=False)
    canvas = FigureCanvas(fig)
    # TODO add a colorbar axis
    gs = gridspec.GridSpec(nfits, nplanes + 1,
                           left=0.08, right=0.91, bottom=0.05, top=0.97,
                           wspace=0.15, hspace=0.2,
                           width_ratios=[1.] * nplanes + [0.1],
                           height_ratios=None)

    cb_ax = fig.add_subplot(gs[:, -1])
    axes = {}
    for i, fit_key in enumerate(fit_labels):
        fit_axes = {k: fig.add_subplot(gs[i, j])
                    for j, k in enumerate(fit_labels)}
        axes[fit_key] = fit_axes

    cube_map = perceptual_rainbow_16.mpl_colormap
    imshow_args = dict(vmax=20, cmap=cube_map)

    for fit_key in usable_fits:
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
                txt = '$\chi^2_\mathrm{{red}}={0:.1f}$'.format(chi_red)
            else:
                txt = '$\chi^2_\mathrm{{red}}={0:.1e}$'.format(chi_red)
            ax.text(0.05, 0.95,
                    txt,
                    ha='left', va='top',
                    backgroundcolor='w',
                    size=9,
                    transform=ax.transAxes)

    for fit_key in fit_labels:
        for plane_key in fit_labels:
            ax = axes[fit_key][plane_key]

            if fit_key == plane_key:
                highlight_color = '#3498db'
                for loc in ['bottom', 'top', 'right', 'left']:
                    ax.spines[loc].set_color(highlight_color)
                    ax.spines[loc].set_linewidth(2.)

            ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(base=1))
            ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(base=0.25))
            ax.xaxis.set_major_locator(xlocators[plane_key])
            ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(base=0.1))

    cb = fig.colorbar(ax=ax, cax=cb_ax, mappable=chi_map)
    cb.set_label(r'$\chi^2$')

    for fit_key in fit_labels.keys()[:-1]:
        for plane_key in fit_labels.keys():
            ax = axes[fit_key][plane_key]
            for tl in ax.get_xmajorticklabels():
                tl.set_visible(False)
            ax.set_xlabel('')

    for fit_key in fit_labels.keys():
        for plane_key in fit_labels.keys()[1:]:
            ax = axes[fit_key][plane_key]
            for tl in ax.get_ymajorticklabels():
                tl.set_visible(False)
            ax.set_ylabel('')

    for fit_key in fit_labels.keys():
        plane_key = fit_labels.keys()[0]
        ax = axes[fit_key][plane_key]
        ax.text(0.00, 1.03, fit_labels[fit_key],
                ha='left', va='baseline',
                transform=ax.transAxes,
                size=10)

    canvas.print_figure(plot_path + ".pdf", format="pdf")


def tabulate_fit_chi(table_path, p, dataset):
    fit_labels = OrderedDict((
        ('lewis', 'Fitting ACS-MS'),
        ('acs_rgb', 'Fitting ACS-RGB'),
        ('acs_all', 'Fitting ACS-ALL'),
        ('oir_all', 'Fitting OIR-ALL'),
        ('ir_rgb', 'Fitting NIR-RGB')))
    nfits = len(fit_labels)
    nplanes = len(fit_labels)

    # p.fit('lewis', ['lewis'], dataset)
    p.fit('acs_rgb', ['acs_rgb'], dataset)
    p.fit('acs_all', ['acs_all'], dataset)
    # p.fit('oir_all', ['oir_all'], dataset)
    p.fit('ir_rgb', ['ir_rgb'], dataset)

    usable_fits = ['acs_rgb', 'acs_all', 'ir_rgb', 'oir_all']
