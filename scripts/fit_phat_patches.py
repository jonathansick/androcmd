#!/usr/bin/env python
# encoding: utf-8
"""
Fit individiual patches in a PHAT footprint.

2015-06-03 - Created by Jonathan Sick
"""


import json


def main():
    args = parse_args()

    # download the star catalog for this brick, if necessary
    download_brick_catalog(args.brick)

    with open(args.json_patch_path, 'r') as f:
        patch_json = json.load(f)

    for patch_num in args.patches:
        # fit this patch
        patch_info = get_patch_info(patch_json, args.brick, patch_num)
        result_hdf5_path = fit_patch(patch_info)

        if args.vodir is not None:
            upload_result(result_hdf5_path, args.vodir)  # FIXME implement


def parse_args():
    pass


def download_brick_catalog(brick):
    """Ensure that the V2 brick star catalog has been downloaded from MAST."""
    pass


def get_patch_info(json, brick, patch_num):
    name = '{0:02d}_{1:03d}'.format(brick, patch_num)
    for d in json:
        if d['patch'] == name:
            return d


def fit_patch(path_info):
    """Fit a patch."""
    pass


def upload_result(result_hdf5_path, vos_dir):
    pass


if __name__ == '__main__':
    main()
