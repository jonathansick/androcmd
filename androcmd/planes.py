#!/usr/bin/env python
# encoding: utf-8
"""
Factories for ColorPlanes used in modelling.

2015-04-22 - Created by Jonathan Sick
"""

from collections import namedtuple, OrderedDict
from starfisher import ColorPlane
from starfisher.pipeline import PlaneBase


Lim = namedtuple('Lim', 'x y')


class BasicPhatPlanes(PlaneBase):
    """Color planes for PHAT data."""
    def __init__(self, **kwargs):
        self._planes = OrderedDict([
            ('f475w_f814w', make_f475w_f814w()),
            ('f475w_f814w_ms', make_f475w_f814w_ms()),
            ('f475w_f814w_rgb', make_f475w_f814w_rgb()),
            ('f475w_f160w', make_f475w_f160w()),
            ('f110w_f160w', make_f110w_f160w()),
        ])
        print "BasicPhatPlanes", kwargs
        super(BasicPhatPlanes, self).__init__(**kwargs)

    @property
    def planes(self):
        return self._planes


class CompletePhatPlanes(PlaneBase):
    """Color planes for PHAT data."""
    def __init__(self, **kwargs):
        self._planes = OrderedDict([
            ('f475w_f814w', make_f475w_f814w()),
            ('f475w_f814w_rgb', make_f475w_f814w_rgb()),
            ('f475w_f110w', make_f475w_f110w()),
            ('f475w_f160w', make_f475w_f160w()),
            ('f814w_f110w', make_f814w_f110w()),
            ('f814w_f160w', make_f814w_f160w()),
            ('f110w_f160w', make_f110w_f160w()),
        ])
        print "CompletePhatPlanes", kwargs
        super(CompletePhatPlanes, self).__init__(**kwargs)

    @property
    def planes(self):
        return self._planes


class RgbPhatPlanes(PlaneBase):
    """Color planes for PHAT data."""
    def __init__(self, **kwargs):
        self._planes = OrderedDict([
            ('f475w_f814w_rgb', make_f475w_f814w_rgb()),
        ])
        print "RgbPhatPlanes", kwargs
        super(RgbPhatPlanes, self).__init__(**kwargs)

    @property
    def planes(self):
        return self._planes


def make_f475w_f814w(dpix=0.05, mag_lim=30.):
    lim = Lim(x=(-1, 5.), y=(25.5, 20.))
    plane = ColorPlane(('F475W', 'F814W'), 'F814W',
                       lim.x,
                       (min(lim.y), max(lim.y)),
                       mag_lim,
                       suffix='f475f814',
                       x_label=r'$\mathrm{F475W}-\mathrm{F814W}$',
                       y_label=r'$\mathrm{F814W}$',
                       dpix=dpix)
    plane.mask_region((3, 5), (28, 25))
    plane.mask_region((3.5, 5), (25, 23))
    plane.mask_region((4, 5), (23, 22.5))
    return plane


def make_f475w_f814w_rgb_hack(dpix=0.05, mag_lim=30.):
    lim = Lim(x=(1.2, 5.), y=(23.5, 20.))
    plane = ColorPlane(('F475W', 'F814W'), 'F814W',
                       lim.x,
                       (min(lim.y), max(lim.y)),
                       mag_lim,
                       suffix='rgbopt',
                       x_label=r'$\mathrm{F475W}-\mathrm{F814W}$',
                       y_label=r'$\mathrm{F814W}$',
                       dpix=dpix,
                       nx=75)  # NOTE auto-calc failed to compute 75
    # plane.mask_region((3, 5), (28, 25))
    # plane.mask_region((3.5, 5), (25, 23))
    # plane.mask_region((4, 5), (23, 22.5))
    return plane


def make_f475w_f814w_rgb(dpix=0.05, mag_lim=30.):
    lim = Lim(x=(1.2, 5.), y=(23.5, 20.))
    plane = ColorPlane(('F475W', 'F814W'), 'F814W',
                       lim.x,
                       (min(lim.y), max(lim.y)),
                       mag_lim,
                       suffix='rgbopt',
                       x_label=r'$\mathrm{F475W}-\mathrm{F814W}$',
                       y_label=r'$\mathrm{F814W}$',
                       dpix=dpix)
    # plane.mask_region((3, 5), (28, 25))
    # plane.mask_region((3.5, 5), (25, 23))
    # plane.mask_region((4, 5), (23, 22.5))
    return plane


def make_f475w_f814w_ms(dpix=0.05, mag_lim=30.):
    lim = Lim(x=(-0.5, 1.), y=(26, 21.))
    plane = ColorPlane(('F475W', 'F814W'), 'F814W',
                       lim.x,
                       (min(lim.y), max(lim.y)),
                       mag_lim,
                       suffix='msopt',
                       x_label=r'$\mathrm{F475W}-\mathrm{F814W}$',
                       y_label=r'$\mathrm{F814W}$',
                       dpix=dpix)
    # plane.mask_region((3, 5), (28, 25))
    # plane.mask_region((3.5, 5), (25, 23))
    # plane.mask_region((4, 5), (23, 22.5))
    return plane


def make_f110w_f160w(dpix=0.05, mag_lim=30.):
    lim = Lim(x=(0.3, 1.3), y=(24., 16.5))
    plane = ColorPlane(('F110W', 'F160W'), 'F160W',
                       lim.x,
                       (min(lim.y), max(lim.y)),
                       mag_lim,
                       suffix='f110f160',
                       x_label=r'$\mathrm{F110W}-\mathrm{F160W}$',
                       y_label=r'$\mathrm{F160W}$',
                       dpix=dpix)
    plane.mask_region((-1., 0.), (22., 16))
    plane.mask_region((0, 0.3), (22., 16))
    plane.mask_region((0.3, 0.7), (20., 16))
    plane.mask_region((0.7, 0.8), (19., 16))
    plane.mask_region((0.8, 0.9), (18., 16))
    plane.mask_region((1.1, 1.5), (28, 21))
    return plane


def make_f475w_f160w(dpix=0.05, mag_lim=30.):
    lim = Lim(x=(-0.8, 8.), y=(25., 17.5))
    plane = ColorPlane(('F475W', 'F160W'), 'F160W',
                       lim.x,
                       (min(lim.y), max(lim.y)),
                       mag_lim,
                       suffix='f475f160',
                       x_label=r'$\mathrm{F475W}-\mathrm{F160W}$',
                       y_label=r'$\mathrm{F110W}$',
                       dpix=dpix)
    return plane


def make_f475w_f110w(dpix=0.05, mag_lim=30.):
    lim = Lim(x=(-0.8, 7.), y=(25., 18.))
    plane = ColorPlane(('F475W', 'F110W'), 'F110W',
                       lim.x,
                       (min(lim.y), max(lim.y)),
                       mag_lim,
                       suffix='f475f110',
                       x_label=r'$\mathrm{F475W}-\mathrm{F110W}$',
                       y_label=r'$\mathrm{F110W}$',
                       dpix=dpix)
    return plane


def make_f814w_f110w(dpix=0.05, mag_lim=30.):
    lim = Lim(x=(-0.1, 1.8), y=(25, 19))
    plane = ColorPlane(('F814W', 'F110W'), 'F110W',
                       lim.x,
                       (min(lim.y), max(lim.y)),
                       mag_lim,
                       suffix='f814f110',
                       x_label=r'$\mathrm{F814W}-\mathrm{F110W}$',
                       y_label=r'$\mathrm{F110W}$',
                       dpix=dpix)
    return plane


def make_f814w_f160w(dpix=0.05, mag_lim=30.):
    lim = Lim(x=(-0.5, 3), y=(24, 17.5))
    plane = ColorPlane(('F814W', 'F160W'), 'F160W',
                       lim.x,
                       (min(lim.y), max(lim.y)),
                       mag_lim,
                       suffix='f814f160',
                       x_label=r'$\mathrm{F814W}-\mathrm{F160W}$',
                       y_label=r'$\mathrm{F160W}$',
                       dpix=dpix)
    return plane
