#!/usr/bin/env python
# encoding: utf-8
"""
Build a JSON file specifying all patches in a PHAT footprint.

2015-05-31 - Created by Jonathan Sick
"""

import argparse
import json

from androcmd.phatpatchfit import build_patches, build_field_patches


def main():
    args = parse_args()

    patches = []

    if args.fields:
        # Patches are PHAT fields
        # This solves problems with filling factors in the brick mosaics
        for brick in args.bricks:
            patches.extend(build_field_patches(brick))
    else:
        # Patches are a given size, segmented from brick images.
        # Problem with this is that the brick images have a very low filling
        # factor
        for brick in args.bricks:
            patches.extend(build_patches(brick, proj_size=args.size))

    with open(args.json_name, 'w') as f:
        json.dump(patches, f, indent=4)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('json_name')
    parser.add_argument('--fields',
                        action='store_true',
                        help='Patches correspond to PHAT fields')
    parser.add_argument('--bricks', nargs='*', type=int,
                        default=range(1, 24))
    parser.add_argument('--size', type=float,
                        help='Projected size of patch in pc on each side',
                        default=100.)
    return parser.parse_args()


if __name__ == '__main__':
    main()
