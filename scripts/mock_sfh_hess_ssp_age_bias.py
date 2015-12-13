"""
Plot mock and residual OIR-ALL hess diagrams for all AST fields, but at a
single age.

2015-12-12 - Created by Jonathan Sick
"""

import os
from collections import namedtuple

import numpy as np
import h5py

import matplotlib as mpl
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.gridspec as gridspec

from palettable.cubehelix import perceptual_rainbow_16
from palettable.colorbrewer.diverging import RdBu_11


HessPlane = namedtuple('HessPlane',
                       'image extent origin y_label x_label x_span y_span')


def main():
    plot('mock_ssp_hess_residuals_5370myr', 'ssp_5370myr_solar')
    plot('mock_ssp_hess_residuals_1380myr', 'ssp_1380myr_solar')


def load_dataset(ast_number, model_name):
    path = os.path.join(os.getenv('STARFISH'),
                        'm{0:d}'.format(ast_number),
                        'm{0:d}.hdf5'.format(ast_number))
    f = h5py.File(path, 'r')

    mock = f['mocksfh'][model_name]['obs_hess']['oir_all']
    print mock.attrs.keys()
    mock_image = np.log10(mock)
    mock_image = np.ma.masked_invalid(mock_image)
    mock_hess = HessPlane(image=mock_image,
                          extent=mock.attrs['extent'],
                          origin=mock.attrs['origin'],
                          y_label=mock.attrs['y_label'],
                          x_label=mock.attrs['x_label'],
                          y_span=mock.attrs['y_span'],
                          x_span=mock.attrs['x_span'])

    diff = f['mocksfh'][model_name]['diff_hess']['oir_all']
    diff_image = np.ma.masked_invalid(diff, copy=True)
    diff_hess = HessPlane(image=diff_image,
                          extent=diff.attrs['extent'],
                          origin=diff.attrs['origin'],
                          y_label=diff.attrs['y_label'],
                          x_label=diff.attrs['x_label'],
                          y_span=diff.attrs['y_span'],
                          x_span=diff.attrs['x_span'])
    f.close()
    return mock_hess, diff_hess


def plot(plot_path, model_name):
    ast_fields = [6, 5, 4, 3]
    d_xticks = 2.
    d_yticks = 1.

    fig = Figure(figsize=(6.5, 4.0), frameon=False)
    canvas = FigureCanvas(fig)
    gs = gridspec.GridSpec(2, len(ast_fields) + 1,
                           left=0.08, right=0.9, bottom=0.15, top=0.90,
                           wspace=0.15, hspace=0.15,
                           width_ratios=[1. for i in ast_fields] + [0.1],
                           height_ratios=None)

    ax_mock_cb = fig.add_subplot(gs[0, -1])
    ax_diff_cb = fig.add_subplot(gs[1, -1])
    axes_mock = []
    axes_diff = []

    for i, ast_field in enumerate(ast_fields):
        ax_mock = fig.add_subplot(gs[0, i])
        axes_mock.append(ax_mock)
        ax_diff = fig.add_subplot(gs[1, i])
        axes_diff.append(ax_diff)

        ax_mock.text(0.5, 1.04,
                     r'AST \#{0:d}'.format(ast_field),
                     va='baseline', ha='center',
                     transform=ax_mock.transAxes)

        mock, diff = load_dataset(ast_field, model_name)

        # Plot the mock image
        im = ax_mock.imshow(mock.image, vmin=0, vmax=1.,
                            cmap=perceptual_rainbow_16.mpl_colormap,
                            extent=mock.extent,
                            origin=mock.origin,
                            norm=None,
                            aspect='auto',
                            interpolation='none')
        if i == 0:
            ax_mock.set_ylabel(mock.y_label)
        else:
            for tl in ax_mock.get_ymajorticklabels():
                tl.set_visible(False)
        for tl in ax_mock.get_xmajorticklabels():
            tl.set_visible(False)
        ax_mock.xaxis.set_major_locator(
            mpl.ticker.MultipleLocator(base=d_xticks))
        ax_mock.yaxis.set_major_locator(
            mpl.ticker.MultipleLocator(base=d_yticks))
        if i == len(ast_fields) - 1:
            cb = fig.colorbar(im, cax=ax_mock_cb, orientation='vertical')
            cb.set_label(r"$\log(N_*)$ Mock")
            # cb.ax.xaxis.set_ticks_position('top')
            # cb.locator = mpl.ticker.MultipleLocator(1.)
            cb.update_ticks()

        # Plot the difference image
        im = ax_diff.imshow(diff.image, vmin=-5, vmax=5,
                            cmap=RdBu_11.mpl_colormap,
                            extent=diff.extent,
                            origin=diff.origin,
                            norm=None,
                            aspect='auto',
                            interpolation='none')
        ax_diff.set_xlabel(mock.x_label)
        if i == 0:
            ax_diff.set_ylabel(mock.y_label)
        else:
            for tl in ax_diff.get_ymajorticklabels():
                tl.set_visible(False)
        ax_diff.xaxis.set_major_locator(
            mpl.ticker.MultipleLocator(base=d_xticks))
        ax_diff.yaxis.set_major_locator(
            mpl.ticker.MultipleLocator(base=d_yticks))
        if i == len(ast_fields) - 1:
            cb = fig.colorbar(im, cax=ax_diff_cb, orientation='vertical')
            cb.set_label(r"$\Delta_\mathrm{mock-model}$ ($N_*$)")
            # cb.ax.xaxis.set_ticks_position('top')
            # cb.locator = mpl.ticker.MultipleLocator(20)
            cb.update_ticks()

    canvas.print_figure(plot_path + ".pdf", format="pdf")


if __name__ == '__main__':
    main()
