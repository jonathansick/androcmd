#!/usr/bin/env python
# encoding: utf-8
"""
Make a CANFAR queue for patch fitting.

>>> from canque import Submission
>>> sub = Submission(user_name, script_path)
>>> sub.add_job('my args', "job.log")
>>> sub.write("jobs.sub")

2015-06-10 - Created by Jonathan Sick
"""

import os
import argparse
import json

import numpy as np
from canque import Submission


def main():
    args = parse_args()

    with open(args.json_patch_path, 'r') as f:
        patch_json = json.load(f)

    patch_numbers = {}
    for brick in args.bricks:
        patch_numbers[brick] = patch_numbers_for_brick(brick, patch_json)

    i = 0
    j = 0
    repeat = True
    while repeat:
        i += 1
        sub = Submission('jonathansick',
                         'androcmd_scripts/patch_fit.sh')
        for brick in args.bricks:
            j += 1
            print len(patch_numbers[brick]),
            selected_patches = select_patches(patch_numbers[brick], args.n)
            print len(patch_numbers[brick])
            add_job(j, sub, brick, selected_patches, args.vodir)

        sub.write(os.path.splitext(args.queue_file)[0]
                  + '_{0:d}.sub'.format(i))

        for brick in args.bricks:
            if len(patch_numbers[brick]) == 0:
                repeat = False


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
    return parser.parse_args()


def patch_numbers_for_brick(brick, patch_json):
    nums = []
    brick_patches = [p for p in patch_json if p['brick'] == brick]

    for patch in brick_patches:
        if patch['brick'] == brick:
            nums.append(int(patch['patch'].split('_')[-1]))
    return nums


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
