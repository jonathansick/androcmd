# encoding: utf-8
"""
Plots comparison maps of ACS-MS and OIR-ALL fitted planes

2015-06-30 - Created by Jonathan Sick
"""

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.gridspec as gridspec

from .galexmap import load_galex_map, setup_galex_axes


def setup_plane_comp_axes():
    basemap = load_galex_map()

    fig = Figure(figsize=(6, 3.5), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(1, 3,
                           left=0.08, right=0.9, bottom=0.1, top=0.95,
                           wspace=0.05, hspace=None,
                           width_ratios=(1, 1, 0.08), height_ratios=None)

    ax_ms = setup_galex_axes(fig, gs[0], basemap)
    ax_oir = setup_galex_axes(fig, gs[1], basemap)
    ax_cb = fig.add_subplot(gs[2])

    for ax in (ax_oir, ax_ms):
        ax.coords[1].set_major_formatter('d.d')
        ax.coords[0].set_major_formatter('hh:mm')
    ax_oir.coords[1].ticklabels.set_visible(False)
    ax_ms.text(0.9, 0.9, 'ACS-MS', transform=ax_ms.transAxes,
               ha='right', va='top',
               zorder=10)
    ax_oir.text(0.9, 0.9, 'OIR-ALL', transform=ax_oir.transAxes,
                ha='right', va='top',
                zorder=10)

    return fig, canvas, ax_ms, ax_oir, ax_cb
