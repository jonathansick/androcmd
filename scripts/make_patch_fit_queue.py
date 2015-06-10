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


import argparse
import json

from canque import Submission


def main():
    args = parse_args()

    with open(args.json_patch_path, 'r') as f:
        patch_json = json.load(f)

    sub = Submission('jonathansick', 'androcmd_scripts/patch_fit.sh')

    job_num = 0
    for brick in args.bricks:
        nums = patch_numbers_for_brick(brick, patch_json)
        while len(nums) > 0:
            create_job(job_num, sub, brick, nums, args.n, args.vodir)
            job_num += 1

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
                        help='Brick number(s)')
    parser.add_argument('--json', dest='json_patch_path',
                        help='Path to patch JSON file')
    parser.add_argument('--vodir',
                        help='VOSpace directory to save results in',
                        default='phat/patches')
    parser.add_argument('--n', type=int,
                        help='Max number of jobs per brick')
    return parser.parse_args()


def patch_numbers_for_brick(brick, patch_json):
    nums = []
    for patch in patch_json:
        if patch['brick'] == brick:
            nums.append(int(patch['patch'].split('_')[-1]))
    return nums


def create_job(job_num, sub, brick, patch_numbers, max_n, vos_dir):
    ns = []
    for i in xrange(max_n):
        if len(patch_numbers) > 0:
            ns.append(str(patch_numbers.pop(0)))
        else:
            break
    job_arg = '{brick:d} {nums} {vos}'.format(
        brick=brick,
        nums=','.join(ns),
        vos=vos_dir)
    sub.add_job(job_arg, "patches_{brick:d}_{job_num:d}.log".format(
        job_num=job_num, brick=brick))


if __name__ == '__main__':
    main()
