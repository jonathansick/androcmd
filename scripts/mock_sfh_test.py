#!/usr/bin/env python
# encoding: utf-8
"""
Run mock star formation history fitting tests.

2015-07-04 - Created by Jonathan Sick
"""

import json
import argparse

import numpy as np

from androcmd.mocksfh.pipeline import RealErrorsThreeZPipeline
from androcmd.mocksfh.pipeline import MockFit


def main():
    args = parse_args()

    with open(args.json_patch_path, 'r') as f:
        patch_json = json.load(f)

    patch_info = get_patch_info(patch_json, args.brick, args.patch)

    isoc = dict(isoc_kind='parsec_CAF09_v1.2S',
                photsys_version='yang')
    kwargs = {}
    kwargs.update(patch_info)
    kwargs['isoc_args'] = isoc
    kwargs['root_dir'] = args.name
    P = PIPELINES[args.pipeline]
    p = P(**kwargs)

    factory = SFH_FACTORIES[args.sfh_name]

    mockfit = MockFit(args.name, factory, p, n_star_amp=True)
    mockfit.make_dataset()
    mockfit.run_fit(args.fit)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('name', type=int,
                        help='Test name')
    parser.add_argument('sfh_name', type=str,
                        choices=SFH_FACTORIES.keys())
    parser.add_argument('--pipeline',
                        default='realistic',
                        choices=PIPELINES.keys())
    parser.add_argument('--brick', type=int, default=23,
                        help='Brick number')
    parser.add_argument('--patch', type=int, default=6,
                        help='Patch number to fit in brick '
                             'can be comma-delimited')
    parser.add_argument('--json', dest='json_patch_path',
                        default='phat_field_patches.json',
                        help='Path to patch JSON file')
    return parser.parse_args()


def ssp_100myr_solar(lockfile):
    ages = 10 ** (lockfile.group_logages - 6.)
    Zs = lockfile.group_metallicities
    n_groups = len(ages)

    # select the single group that matches the description
    i = np.argmin(np.hypot(ages - 100., Zs - 0.019))

    sfhs = np.zeros(n_groups, dtype=np.float)
    sfhs[i] = 10000
    return sfhs


SFH_FACTORIES = {
    'ssp_100myr_solar': ssp_100myr_solar,
}


PIPELINES = {
    'realistic': RealErrorsThreeZPipeline,
}


def get_patch_info(json, brick, patch_num):
    name = '{0:02d}_{1:03d}'.format(brick, patch_num)
    for d in json:
        if d['patch'] == name:
            return d


if __name__ == '__main__':
    main()
