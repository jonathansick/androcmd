#!/usr/bin/env python
# encoding: utf-8
"""
Plotting utilities for androcmd
"""

import matplotlib as mpl
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
