#!/usr/bin/env python
# encoding: utf-8
"""
Plots of cumulative mass formation for mock trials with tau star formation
histories.

2016-01-18 - Created by Jonathan Sick
"""
import os
import h5py
import numpy as np
from astropy.table import Table
from starfisher.sfh import marginalize_sfh_metallicity
from palettable.colorbrewer.qualitative import Dark2_6

import matplotlib as mpl
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.gridspec as gridspec


def main():
    plot('mock_sfh_tau_cumulative_mass', show_logage=False)
    plot('mock_sfh_tau_cumulative_mass_logage', show_logage=True)


def plot(plot_path, show_logage=False):
    # Star formation histories
    model_list = ['tau_0.1_solar',
                  'tau_0.5_solar',
                  'tau_1.0_solar',
                  'tau_5.0_solar',
                  'tau_10.0_solar',
                  'tau_20.0_solar',
                  'tau_50.0_solar',
                  'tau_100.0_solar']
    model_taus = np.array([0.1, 0.5, 1.0, 5.0, 10.0, 20.0, 50.0, 100.0])

    # Fitting experiements (AST fields, and errorless/ deep)
    root_path = os.getenv('STARFISH')
    all_colors = Dark2_6.mpl_colors
    realistic_linestyle = {'ls': '-', 'lw': 1}
    errorless_linestyle = {'ls': '--', 'dashes': (5, 2), 'lw': 1}
    # errorless_linestyle = {'ls': '--', 'lw': 1}
    experiments = [
        [os.path.join(root_path, 'm3', 'm3.hdf5'),
         'oir_all',
         r'AST \#3',
         all_colors[0],
         realistic_linestyle],
        [os.path.join(root_path, 'm4', 'm4.hdf5'),
         'oir_all',
         r'AST \#4',
         all_colors[1],
         realistic_linestyle],
        [os.path.join(root_path, 'm5', 'm5.hdf5'),
         'oir_all',
         r'AST \#5',
         all_colors[2],
         realistic_linestyle],
        [os.path.join(root_path, 'm6', 'm6.hdf5'),
         'oir_all',
         r'AST \#6',
         all_colors[3],
         realistic_linestyle],
        [os.path.join(root_path, 'idealall', 'idealall.hdf5'),
         'oir_all',
         r'Errorless',
         all_colors[4],
         errorless_linestyle],
        [os.path.join(root_path, 'idealall', 'idealall.hdf5'),
         'oir_all_28',
         r'Deep Errorless',
         all_colors[5],
         errorless_linestyle],
    ]

    nx = 4
    ny = 2
    assert nx * ny == len(model_list)

    fig = Figure(figsize=(6.5, 4.5), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(ny, nx,
                           left=0.09, right=0.8, bottom=0.12, top=0.95,
                           wspace=0.15, hspace=None,
                           width_ratios=None, height_ratios=None)
    for i, (model_name, model_tau) in enumerate(zip(model_list, model_taus)):
        iy = i / nx
        ix = i % nx
        ax = fig.add_subplot(gs[iy, ix])
        ax.text(0.5, 1.02, r'$\tau = {0:.1f}$'.format(model_tau),
                ha='center', va='bottom', transform=ax.transAxes)

        for e in experiments:
            exp_path, fit_key, exp_label, c, ls = e
            logage, model_cmass, fit_cmass = extract_cumulative_mass_function(
                exp_path, model_name, fit_key)

            if show_logage:
                A = logage
            else:
                A = 10. ** (logage - 9.)
            ax.plot(A,
                    fit_cmass,
                    label=exp_label,
                    c=c,
                    **ls)
            ax.set_ylim(-0.05, 1.05)
            if show_logage:
                ax.set_xlim(6.5, 10.2)
            else:
                ax.set_xlim(0., 12.)
        ax.plot(A,
                model_cmass,
                ls='-',
                c='0.7',
                alpha=1.,
                lw=4,
                zorder=-1,
                label='Model')
        if iy == 0 and ix == 3:
            ax.legend(bbox_to_anchor=(1.05, 1), loc=2, fontsize=7,
                      frameon=False)

        if ix > 0:
            for tl in ax.get_ymajorticklabels():
                tl.set_visible(False)
        else:
            ax.set_ylabel(r'$M(t_\mathrm{L}>A) / \sum M$')
        if iy < ny - 1:
            for tl in ax.get_xmajorticklabels():
                tl.set_visible(False)
        else:
            if show_logage:
                ax.set_xlabel(r'$\log_{10} (A~\mathrm{Gyr}^{-1})$')
                ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(base=1))
            else:
                ax.set_xlabel(r'$A$ (Gyr)')
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(8)

    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    canvas.print_figure(plot_path + ".pdf", format="pdf")


def extract_cumulative_mass_function(dataset_path, model_key, fit_key):
    """Build a cumulative mass distribution function for a single plane,
    for a single of tau-model fit.
    """
    def _compute_cmass(d):
        d = np.array(d)
        loga = d['log(age)']
        cmass = np.cumsum(d['mass'][::-1])[::-1]
        cmass = cmass / cmass.max()  # normalize
        return loga, cmass

    def _marginalize_sfh(d):
        """Marginalize the SFH table"""
        t = Table(np.array(d))
        marginalized_t = marginalize_sfh_metallicity(t)
        return marginalized_t

    f = h5py.File(dataset_path)
    model_sfh = f['mocksfh'][model_key]['mock_sfh_marginal']
    fit_sfh = _marginalize_sfh(f['mocksfh'][model_key]['sfh'][fit_key])
    logage, model_cmass = _compute_cmass(model_sfh)
    logage, fit_cmass = _compute_cmass(fit_sfh)
    f.close()

    return logage, model_cmass, fit_cmass


if __name__ == '__main__':
    main()
