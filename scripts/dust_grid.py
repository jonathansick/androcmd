#!/usr/bin/env python
# encoding: utf-8
"""
Make a grid of synths for a set of attenuations.

2015-04-30 - Created by Jonathan Sick
"""

import argparse

import numpy as np

from starfisher.pipeline import PipelineBase
from androcmd.planes import BasicPhatPlanes
from androcmd.phatpipeline import (
    SolarZIsocs, SolarLockfile,
    PhatGaussianDust, PhatCrowding)
# from androcmd.phatpipeline import PhatCatalog


def main():
    args = parse_args()

    av_grid = np.arange(0., args.max_av, args.delta_av)
    for av in av_grid:
        run_pipeline(brick=args.brick, av=av)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Grid of synths for a set of Av")
    parser.add_argument('brick', type=int)
    parser.add_argument('max_av', type=float, default=1.5)
    parser.add_argument('delta_av', type=float, default=0.1)
    return parser.parse_args()


def run_pipeline(brick=23, av=0., run_fit=False):
    # dataset = PhatCatalog(brick)
    pipeline = Pipeline(root_dir="b{0:d}_{1:.2f}".format(brick, av),
                        young_av=av, old_av=av, av_sigma_ratio=0.25,
                        isoc_args=dict(isoc_kind='parsec_CAF09_v1.2S',
                                       photsys_version='yang'))
    print(pipeline)
    print('av {0:.1f} done'.format(av))


class Pipeline(BasicPhatPlanes, SolarZIsocs,
               SolarLockfile, PhatGaussianDust, PhatCrowding, PipelineBase):
    """A pipeline for fitting PHAT bricks with solar metallicity isochrones."""
    def __init__(self, **kwargs):
        super(Pipeline, self).__init__(**kwargs)


if __name__ == '__main__':
    main()
