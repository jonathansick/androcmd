#!/usr/bin/env python
# encoding: utf-8
"""
Build a JSON file specifying all patches in a PHAT footprint.

2015-05-31 - Created by Jonathan Sick
"""

import argparse
import json

from androcmd.phatpatchfit import build_patches


def main():
    args = parse_args()

    patches = []

    for brick in args.bricks:
        patches.extend(build_patches(brick, proj_size=args.size))

    with open(args.json_name, 'w') as f:
        json.dump(patches, f, indent=4)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('json_name')
    parser.add_argument('--bricks', nargs='*', type=int,
                        default=range(2, 24))
    parser.add_argument('--size', type=float,
                        help='Projected size of patch in pc on each side',
                        default=100.)
    return parser.parse_args()


if __name__ == '__main__':
    main()
