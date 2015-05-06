#!/usr/bin/env python
# encoding: utf-8
"""
Modelling a brick's CMD with different choices of circumstellar dust.

2015-05-06 - Created by Jonathan Sick
"""

import argparse

from starfisher.pipeline import PipelineBase
from androcmd.planes import BasicPhatPlanes
from androcmd.phatpipeline import (
    SolarZIsocs, SolarLockfile,
    NoDust, PhatCrowding)
from androcmd.phatpipeline import PhatCatalog


def main():
    args = parse_args()

    pipeline = make_pipeline(args.brick, args.dustm, args.dustc)

    fit_planes = {'f475w_f160w': ['f475w_f160w'],
                  'rgb': ['f475w_f814w_rgb'],
                  'ms': ['f475w_f814w_ms']}
    if len(fit_planes) > 0:
        dataset = PhatCatalog(args.brick)
    for fit_name in args.fit:
        plane = fit_planes[fit_name]
        pipeline.fit(fit_name, plane, dataset)


def make_pipeline(brick, dustM, dustC):
    isoc_args = dict(isoc_kind='parsec_CAF09_v1.2S',
                     photsys_version='yang')
    isoc_args['dust_sourceM'] = dustM
    isoc_args['dust_sourceC'] = dustC

    pipeline = Pipeline(root_dir="b{0:d}_{1}_{2}".format(brick,
                                                         dustM,
                                                         dustC),
                        isoc_args=isoc_args)

    return pipeline


def parse_args():
    parser = argparse.ArgumentParser(
        description="Fitting with different circumstellar dust models")
    parser.add_argument('brick', type=int)
    parser.add_argument('--fit', nargs='*', default=list())
    parser.add_argument('--dustm', choices=['nodustM', 'sil', 'AlOx',
                                            'dpmod60alox40', 'dpmod'],
                        default='nodustM')
    parser.add_argument('--dustc', choices=['nodustC', 'gra', 'AMC',
                                            'AMCSIC15'],
                        default='gra')
    return parser.parse_args()


class Pipeline(BasicPhatPlanes, SolarZIsocs,
               SolarLockfile, NoDust, PhatCrowding, PipelineBase):
    """A pipeline for fitting PHAT bricks with solar metallicity isochrones."""
    def __init__(self, **kwargs):
        super(Pipeline, self).__init__(**kwargs)


if __name__ == '__main__':
    main()
