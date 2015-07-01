#!/usr/bin/env python
# encoding: utf-8
"""
Make a CANFAR queue for patch fitting.

2015-06-10 - Created by Jonathan Sick
"""

import argparse
import json
import h5py

import numpy as np
from canque import Submission


def main():
    args = parse_args()

    with open(args.json_patch_path, 'r') as f:
        patch_json = json.load(f)

    if args.existing is not None:
        dataset = h5py.File(args.existing, 'r')
    else:
        dataset = None

    patch_numbers = {}
    for brick in args.bricks:
        patch_numbers[brick] = patch_numbers_for_brick(brick, patch_json,
                                                       dataset)

    i = 0
    sub = Submission('jonathansick',
                     'androcmd_scripts/patch_fit.sh')
    for brick in args.bricks:
        while len(patch_numbers[brick]) > 0:
            selected_patches = select_patches(patch_numbers[brick], args.n)
            add_job(i, sub, brick, selected_patches, args.vodir)
            i += 1

    sub.write(args.queue_file)


def parse_args():
    parser = argparse.ArgumentParser(
        description='e.g.:\n\nmake_patch_fit_queue.py brick_23_queue.sub '
                    '--bricks 23 '
                    '--json ~/Downloads/patches.json --n 30')
    parser.add_argument('queue_file',
                        help='Output path of queue submission file')
    parser.add_argument('--bricks', type=int,
                        nargs='*',
                        default=range(2, 24),
                        help='Brick number(s)')
    parser.add_argument('--json', dest='json_patch_path',
                        help='Path to patch JSON file')
    parser.add_argument('--vodir',
                        help='VOSpace directory to save results in',
                        default='phat/patches')
    parser.add_argument('--n', type=int,
                        help='Max number of patches per jobs')
    parser.add_argument('--existing',
                        help='HDF5 results file with existing patches to skip',
                        default=None)
    return parser.parse_args()


def patch_numbers_for_brick(brick, patch_json, dataset):
    nums = []
    brick_patches = [p for p in patch_json if p['brick'] == brick]

    for patch in brick_patches:
        if patch['brick'] == brick:
            nums.append(int(patch['patch'].split('_')[-1]))

    # remove already-computed patches
    if dataset is not None:
        existing = set(existing_patches_for_brick(brick, dataset))
        nums = list(set(nums) - existing)

    return nums


def existing_patches_for_brick(brick, dataset):
    patch_numbers = []
    patches = [k.split('_') for k in dataset['patches'].keys()]
    for brick_str, patch_str in patches:
        if int(brick_str) == brick:
            patch_numbers.append(int(patch_str))
    return patch_numbers


def select_patches(patch_numbers, n_choose):
    selected_nums = []
    if n_choose > len(patch_numbers):
        n_choose = len(patch_numbers)
    for i in xrange(n_choose):
        rand_int = np.random.randint(0,
                                     high=len(patch_numbers),
                                     size=(1,))[0]
        selected_nums.append(patch_numbers.pop(rand_int))
    return [str(i) for i in selected_nums]


def add_job(job_num, sub, brick, patch_numbers, vos_dir):
    job_arg = '{brick:d} {nums} {vos}'.format(
        brick=brick,
        nums=','.join(patch_numbers),
        vos=vos_dir)
    sub.add_job(job_arg, "patches_{brick:d}_{job_num:d}.log".format(
        job_num=job_num, brick=brick))


if __name__ == '__main__':
    main()
