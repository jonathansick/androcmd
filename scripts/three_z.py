#!/usr/bin/env python
# encoding: utf-8
"""
Modelling with three metallicity groups.
"""

import argparse

from starfisher.pipeline import PipelineBase
from androcmd.planes import BasicPhatPlanes
from androcmd.phatpipeline import (
    ExtendedSolarIsocs, ExtendedSolarLockfile,
    NoDust, PhatCrowding)
from androcmd.phatpipeline import PhatCatalog


def main():
    args = parse_args()

    pipeline = Pipeline(root_dir="b{0:d}_threez".format(args.brick),
                        isoc_args=dict(isoc_kind='parsec_CAF09_v1.2S',
                                       photsys_version='yang'))
    fit_planes = {'f475w_f160w': ['f475w_f160w'],
                  'rgb': ['f475w_f814w_rgb'],
                  'ms': ['f475w_f814w_ms']}
    if len(fit_planes) > 0:
        dataset = PhatCatalog(args.brick)
    for fit_name in args.fit:
        plane = fit_planes[fit_name]
        pipeline.fit(fit_name, plane, dataset)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Multi-metallicity fitting")
    parser.add_argument('brick', type=int)
    parser.add_argument('--fit', nargs='*', default=list())
    return parser.parse_args()


class Pipeline(BasicPhatPlanes,
               ExtendedSolarIsocs, ExtendedSolarLockfile,
               NoDust, PhatCrowding, PipelineBase):
    """A pipeline for fitting PHAT bricks with solar metallicity isochrones."""
    def __init__(self, **kwargs):
        print "MultiZPhatPipeline", kwargs
        super(Pipeline, self).__init__(**kwargs)


if __name__ == '__main__':
    main()
