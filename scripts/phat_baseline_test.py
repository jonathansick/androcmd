#!/usr/bin/env python
# encoding: utf-8
"""
Grid computation of dust attenuation for old vs. young stellar populations.

2015-05-12 - Created by Jonathan Sick
"""

import argparse


from androcmd.phatpipeline import PhatCatalog
from androcmd.baselineexp import SolarZPipeline, ThreeZPipeline


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

    if args.plot_hess is not None:
        from androcmd.baselineexp import plot_fit_hess_grid
        dataset = PhatCatalog(args.brick)
        plot_fit_hess_grid(args.plot_hess, pipeline, dataset)

    if args.plot_diff is not None:
        from androcmd.baselineexp import plot_diff_hess_grid
        dataset = PhatCatalog(args.brick)
        plot_diff_hess_grid(args.plot_diff, pipeline, dataset)

    if args.chi_table is not None:
        from androcmd.baselineexp import tabulate_fit_chi
        dataset = PhatCatalog(args.brick)
        tabulate_fit_chi(args.chi_table, pipeline, dataset)


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
    parser.add_argument('--plot-hess', default=None)
    parser.add_argument('--plot-diff', default=None)
    parser.add_argument('--chi-table', default=None)
    return parser.parse_args()


if __name__ == '__main__':
    main()
