#!/usr/bin/env python
# encoding: utf-8
"""
Run mock star formation history fitting tests.

2015-07-04 - Created by Jonathan Sick
"""

import json
import argparse
from functools import partial

import numpy as np

from androcmd.mocksfh.pipeline import (RealErrorsThreeZPipeline,
                                       IdealizedThreeZPipeline)
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
    kwargs['synth_config_only'] = True
    kwargs['synth_config_path'] = 'testpop/synth.txt'
    P = PIPELINES[args.pipeline]
    p = P(**kwargs)

    factory = SFH_FACTORIES[args.sfh_name]

    mockfit = MockFit(args.name, factory, p, n_star_amp=True)
    mockfit.make_dataset()
    mockfit.run_fit(args.fit, n_synth_cpu=args.n_synth_cpu)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('name',
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
    parser.add_argument('--ncpu', type=int, default=4,
                        dest='n_synth_cpu',
                        help='Number of CPUs to run synth with')
    parser.add_argument('--fit', nargs='*',
                        help='Names of fitted planes')
    return parser.parse_args()


def ssp_solar(lockfile, age_myr=100):
    ages = 10 ** (lockfile.group_logages - 6.)
    Zs = lockfile.group_metallicities
    n_groups = len(ages)

    # select the single group that matches the description
    i = np.argmin(np.hypot(ages - age_myr, Zs - 0.019))

    sfhs = np.zeros(n_groups, dtype=np.float)
    sfhs[i] = 1.
    return sfhs


SFH_FACTORIES = dict()
for myr in (100, 250, 500):
    key = 'ssp_{0:d}myr_solar'.format(myr)
    SFH_FACTORIES[key] = partial(ssp_solar, age_myr=myr)
for gyr in range(1, 13):
    key = 'ssp_{0:d}gyr_solar'.format(gyr)
    SFH_FACTORIES[key] = partial(ssp_solar, age_myr=gyr * 1e3)


PIPELINES = {
    'realistic': RealErrorsThreeZPipeline,
    'ideal': IdealizedThreeZPipeline,
}


def get_patch_info(json, brick, patch_num):
    name = '{0:02d}_{1:03d}'.format(brick, patch_num)
    for d in json:
        if d['patch'] == name:
            return d


if __name__ == '__main__':
    main()