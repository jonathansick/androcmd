#!/usr/bin/env python
# encoding: utf-8
"""
Fit all of M31, divided into patches.

Relies on having a JSON file specifying each field. Spec is:

    {patch: <int> patch number,
     poly: <list of (ra,dec) tuples> footprint,
     brick: <int> phat brick number that photometry belongs to,
     ra0: central RA of the patch,
     dec0: central DEC of the patch,
     area: deprojected area of patch in pc^2
     area_proj: projected area of patch in arcsec^2
    }

2015-05-28 - Created by Jonathan Sick
"""

import argparse
import json

# FIXME implement these
from androcmd.phatpatchfit import PatchCatalog
# from androcmd.phatpatchfit import SolarZPipeline, ThreeZPipeline
from androcmd.phatpatchfit import ThreeZPipeline


def main():
    args = parse_args()

    with open(args.patch_file_path) as f:
        patches = json.load(f)

    model_name = "{0:06d}".format(args.patch_number)
    patch = patches[args.patch_number]

    # if args.pipeline == 'solarz':
    #     # Use the single-Z solar pipeline
    #     # Pipeline = SolarZPipeline
    #     pass
    # elif args.pipeline == 'threez':
    #     # Use the three-metallicity track pipeline
    #     # Pipeline = ThreeZPipeline
    #     pass
    Pipeline = ThreeZPipeline
    isoc = dict(isoc_kind='parsec_CAF09_v1.2S',
                photsys_version='yang')
    pipeline = Pipeline(root_dir=model_name,
                        isoc_args=isoc,
                        **patch)
    dataset = PatchCatalog(**patch)

    for fit in args.fit:
        pipeline.fit(fit, [fit], dataset)

    # FIXME perhaps we also want to pre-reduce the results into an HDF5 file
    # for uploading/plotting


def parse_args(arg1):
    parser = argparse.ArgumentParser(
        description="Model a brick with differential old/young dust.")
    parser.add_argument('patch_file_path')
    parser.add_argument('patch_number', type=int)
    parser.add_argument('--fit',
                        choices=['lewis', 'acs_rgb', 'acs_all',
                                 'oir_all', 'ir_rgb'],
                        default=['oir_all'])
    parser.add_argument('--pipeline',
                        choices=['solarz', 'threez'],
                        default='threez')
    return parser.parse_args()


if __name__ == '__main__':
    main()
