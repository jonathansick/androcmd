#!/usr/bin/env python
# encoding: utf-8
"""
Identify which PHAT brick+field each AST field corresponds to.

2015-07-08 - Created by Jonathan Sick
"""

import numpy as np
from sklearn.cluster import KMeans
# from astroML.stats import binned_statistic

# import matplotlib as mpl
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.gridspec as gridspec
from matplotlib.patches import Polygon
from matplotlib.path import Path

from m31hst.phatast import load_phat_ast_table

from androcmd.phatpatchfit.galexmap import (load_galex_map, setup_galex_axes,
                                            plot_patch_footprints)
from androcmd.phatpatchfit.pipeline import load_field_patches


def main():
    fields = load_field_patches()
    ast_centers = load_ast_centers()
    matches = match_fields(fields, ast_centers)
    print matches
    plot_ast_fields(fields, matches, ast_centers=ast_centers)


def load_ast_centers():
    t = load_phat_ast_table()
    km = KMeans(n_clusters=6)
    xy = np.vstack((t['ra'], t['dec'])).T
    km.fit(xy)
    centers = km.cluster_centers_

    srt = np.argsort(centers[:, 1])
    return centers[srt, :]


def plot_ast_fields(fields, matches, ast_centers=None):
    fig = Figure(figsize=(3.5, 3.5), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(1, 1,
                           left=0.15, right=0.95, bottom=0.15, top=0.95,
                           wspace=None, hspace=None,
                           width_ratios=None, height_ratios=None)
    basemap = load_galex_map()
    ax = setup_galex_axes(fig, gs[0], basemap)
    plot_patch_footprints(ax, alpha=0.8, edgecolor='dodgerblue')
    for n, m in matches.iteritems():
        footprint = np.array(m['poly'])
        patch = Polygon(footprint, closed=True,
                        transform=ax.get_transform('world'),
                        facecolor='y', alpha=1,
                        edgecolor='k', lw=0.5, zorder=10)
        ax.add_patch(patch)
        x = footprint[:, 0].mean()
        y = footprint[:, 1].mean()
        ax.annotate('{0:d}'.format(n), xy=(x, y),
                    xycoords=ax.get_transform('world'),
                    xytext=(3, -3), textcoords="offset points",
                    size=8,
                    bbox=dict(boxstyle="round",
                              fc=(1., 1., 1., 0.8),
                              edgecolor='None'))
    if ast_centers is not None:
        ax.scatter(ast_centers[:, 0], ast_centers[:, 1],
                   marker='*', c='y',
                   transform=ax.get_transform('world'))
    gs.tight_layout(fig, pad=1.08, h_pad=None, w_pad=None, rect=None)
    ax.coords[0].ticklabels.set_size(11)
    ax.coords[1].ticklabels.set_size(11)
    canvas.print_figure("phat_ast_fields.pdf", format="pdf")


def match_fields(fields, centers):
    matches = {}
    for i in xrange(centers.shape[0]):
        for field in fields:
            # poly = Polygon(field['poly'])
            poly = Path(field['poly'])
            if poly.contains_point((centers[i, 0], centers[i, 1])):
                match = {'brick': field['brick'],
                         'field': field['field'],
                         'ra': centers[i, 0],
                         'dec': centers[i, 1],
                         'poly': field['poly']}
                matches[i + 1] = match
    return matches


if __name__ == '__main__':
    main()
