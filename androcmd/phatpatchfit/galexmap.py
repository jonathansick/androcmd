# encoding: utf-8
"""
Plotting services with galex basemaps.

2015-06-30 - Created by Jonathan Sick
"""

from collections import namedtuple

import numpy as np
import matplotlib as mpl
# from matplotlib.figure import Figure
# from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
# import matplotlib.gridspec as gridspec
import wcsaxes
import astropy.io.fits

from .pipeline import load_field_footprints

BaseMap = namedtuple('BaseMap', 'image wcs xlim ylim vmin vmax')


def load_galex_map(ref_path='h_m31-nd-int.fits'):
    with astropy.io.fits.open(ref_path) as f:
        header = f[0].header
        base_image = f[0].data
        wcs = wcsaxes.WCS(header)
    basemap = BaseMap(image=np.log10(base_image),
                      wcs=wcs,
                      xlim=(500, 3000),
                      ylim=(2800, 6189),
                      vmin=-2,
                      vmax=-1)
    return basemap


def setup_galex_axes(fig, gs_span, basemap):
    ax = fig.add_subplot(gs_span, projection=basemap.wcs)
    ax.set_xlim(-0.5, basemap.image.shape[1] - 0.5)
    ax.set_ylim(-0.5, basemap.image.shape[0] - 0.5)
    ax.imshow(basemap.image,
              cmap=mpl.cm.gray_r, vmin=basemap.vmin, vmax=basemap.vmax,
              zorder=-10,
              origin='lower')
    ax.set_xlim(500, 3000)
    ax.set_ylim(2800, 6189)
    ax.coords[1].set_major_formatter('d.d')
    ax.coords[0].set_major_formatter('hh:mm')
    ax.coords[0].set_separator((r'$^h$', "'", '"'))
    ax.coords[0].ticklabels.set_size(8)
    plot_patch_footprints(ax)
    return ax


def plot_patch_footprints(ax):
    # Plot phat footprints
    for footprint in load_field_footprints():
        patch = mpl.patches.Polygon(footprint, closed=True,
                                    transform=ax.get_transform('world'),
                                    facecolor='None', alpha=0.1,
                                    edgecolor='k', lw=0.5)
        ax.add_patch(patch)
