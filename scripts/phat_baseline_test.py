#!/usr/bin/env python
# encoding: utf-8
"""
Grid computation of dust attenuation for old vs. young stellar populations.

2015-05-12 - Created by Jonathan Sick
"""

import argparse

from starfisher.pipeline import PipelineBase
from androcmd.planes import BaselineTestPhatPlanes
from androcmd.phatpipeline import (PhatCatalog, SolarZIsocs, SolarLockfile,
                                   LewisBrickDust, PhatCrowding,
                                   ExtendedSolarIsocs, ExtendedSolarLockfile)


def main():
    args = parse_args()

    if args.pipeline == 'solarz':
        # Use the single-Z solar pipeline
        Pipeline = SolarZPipeline
    elif args.pipeline == 'threez':
        # Use the three-metallicity track pipeline
        Pipeline = ThreeZPipeline

    isoc = dict(isoc_kind='parsec_CAF09_v1.2S',
                photsys_version='yang')
    pipeline = Pipeline(brick=23,
                        root_dir=args.model_name,
                        isoc_args=isoc)

    if args.fit is not None:
        dataset = PhatCatalog(args.brick)
        pipeline.fit(args.fit, [args.fit], dataset)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Model a brick with differential old/young dust.")
    parser.add_argument('model_name')
    parser.add_argument('brick', type=int)
    parser.add_argument('--fit',
                        choices=['lewis', 'acs_rgb', 'acs_all',
                                 'oir_all', 'ir_rgb'],
                        default=None)
    parser.add_argument('--pipeline',
                        choices=['solarz', 'threez'],
                        default='solarz')
    return parser.parse_args()


class SolarZPipeline(BaselineTestPhatPlanes, SolarZIsocs,
                     SolarLockfile, LewisBrickDust, PhatCrowding,
                     PipelineBase):
    """Pipeline for the baseline test with only solar metallicity."""
    def __init__(self, **kwargs):
        super(SolarZPipeline, self).__init__(**kwargs)


class ThreeZPipeline(BaselineTestPhatPlanes, ExtendedSolarIsocs,
                     ExtendedSolarLockfile, LewisBrickDust, PhatCrowding,
                     PipelineBase):
    """Pipeline for baseline test with three metallicity tracks."""
    def __init__(self, **kwargs):
        super(ThreeZPipeline, self).__init__(**kwargs)


if __name__ == '__main__':
    main()
