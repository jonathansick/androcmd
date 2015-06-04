#!/usr/bin/env python
# encoding: utf-8
"""
Fit individiual patches in a PHAT footprint.

2015-06-03 - Created by Jonathan Sick
"""

import os
import argparse
import json
import subprocess

from androcmd.phatpatchfit import PatchCatalog, ThreeZPipeline


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
    parser = argparse.ArgumentParser()
    parser.add_argument('brick', type=int,
                        help='Brick number')
    parser.add_argument('--patches', type=int, nargs='*',
                        help='Patch number(s) to fit in brick')
    parser.add_argument('--json', dest='json_patch_path',
                        help='Path to patch JSON file')
    parser.add_argument('--vodir',
                        help='VOSpace directory to save results in')
    return parser.parse_arguments()


def download_brick_catalog(brick):
    """Ensure that the V2 brick star catalog has been downloaded from MAST."""
    urls = {1: 'http://archive.stsci.edu/pub/hlsp/phat/brick01/hlsp_phat_hst_wfc3-uvis-acs-wfc-wfc3-ir_12058-m31-b01_f275w-f336w-f475w-f814w-f110w-f160w_v2_st.fits',  # NOQA
            2: 'http://archive.stsci.edu/pub/hlsp/phat/brick02/hlsp_phat_hst_wfc3-uvis-acs-wfc-wfc3-ir_12073-m31-b02_f275w-f336w-f475w-f814w-f110w-f160w_v2_st.fits',  # NOQA
            3: 'http://archive.stsci.edu/pub/hlsp/phat/brick03/hlsp_phat_hst_wfc3-uvis-acs-wfc-wfc3-ir_12109-m31-b03_f275w-f336w-f475w-f814w-f110w-f160w_v2_st.fits',  # NOQA
            4: 'http://archive.stsci.edu/pub/hlsp/phat/brick04/hlsp_phat_hst_wfc3-uvis-acs-wfc-wfc3-ir_12107-m31-b04_f275w-f336w-f475w-f814w-f110w-f160w_v2_st.fits',  # NOQA
            5: 'http://archive.stsci.edu/pub/hlsp/phat/brick05/hlsp_phat_hst_wfc3-uvis-acs-wfc-wfc3-ir_12074-m31-b05_f275w-f336w-f475w-f814w-f110w-f160w_v2_st.fits',  # NOQA
            6: 'http://archive.stsci.edu/pub/hlsp/phat/brick06/hlsp_phat_hst_wfc3-uvis-acs-wfc-wfc3-ir_12105-m31-b06_f275w-f336w-f475w-f814w-f110w-f160w_v2_st.fits',  # NOQA
            7: 'http://archive.stsci.edu/pub/hlsp/phat/brick07/hlsp_phat_hst_wfc3-uvis-acs-wfc-wfc3-ir_12113-m31-b07_f275w-f336w-f475w-f814w-f110w-f160w_v2_st.fits',  # NOQA
            8: 'http://archive.stsci.edu/pub/hlsp/phat/brick08/hlsp_phat_hst_wfc3-uvis-acs-wfc-wfc3-ir_12075-m31-b08_f275w-f336w-f475w-f814w-f110w-f160w_v2_st.fits',  # NOQA
            9: 'http://archive.stsci.edu/pub/hlsp/phat/brick09/hlsp_phat_hst_wfc3-uvis-acs-wfc-wfc3-ir_12057-m31-b09_f275w-f336w-f475w-f814w-f110w-f160w_v2_st.fits',  # NOQA
            10: 'http://archive.stsci.edu/pub/hlsp/phat/brick10/hlsp_phat_hst_wfc3-uvis-acs-wfc-wfc3-ir_12111-m31-b10_f275w-f336w-f475w-f814w-f110w-f160w_v2_st.fits',  # NOQA
            11: 'http://archive.stsci.edu/pub/hlsp/phat/brick11/hlsp_phat_hst_wfc3-uvis-acs-wfc-wfc3-ir_12115-m31-b11_f275w-f336w-f475w-f814w-f110w-f160w_v2_st.fits',  # NOQA
            12: 'http://archive.stsci.edu/pub/hlsp/phat/brick12/hlsp_phat_hst_wfc3-uvis-acs-wfc-wfc3-ir_12071-m31-b12_f275w-f336w-f475w-f814w-f110w-f160w_v2_st.fits',  # NOQA
            13: 'http://archive.stsci.edu/pub/hlsp/phat/brick13/hlsp_phat_hst_wfc3-uvis-acs-wfc-wfc3-ir_12114-m31-b13_f275w-f336w-f475w-f814w-f110w-f160w_v2_st.fits',  # NOQA
            14: 'http://archive.stsci.edu/pub/hlsp/phat/brick14/hlsp_phat_hst_wfc3-uvis-acs-wfc-wfc3-ir_12072-m31-b14_f275w-f336w-f475w-f814w-f110w-f160w_v2_st.fits',  # NOQA
            15: 'http://archive.stsci.edu/pub/hlsp/phat/brick15/hlsp_phat_hst_wfc3-uvis-acs-wfc-wfc3-ir_12056-m31-b15_f275w-f336w-f475w-f814w-f110w-f160w_v2_st.fits',  # NOQA
            16: 'http://archive.stsci.edu/pub/hlsp/phat/brick16/hlsp_phat_hst_wfc3-uvis-acs-wfc-wfc3-ir_12106-m31-b16_f275w-f336w-f475w-f814w-f110w-f160w_v2_st.fits',  # NOQA
            17: 'http://archive.stsci.edu/pub/hlsp/phat/brick17/hlsp_phat_hst_wfc3-uvis-acs-wfc-wfc3-ir_12059-m31-b17_f275w-f336w-f475w-f814w-f110w-f160w_v2_st.fits',  # NOQA
            18: 'http://archive.stsci.edu/pub/hlsp/phat/brick18/hlsp_phat_hst_wfc3-uvis-acs-wfc-wfc3-ir_12108-m31-b18_f275w-f336w-f475w-f814w-f110w-f160w_v2_st.fits',  # NOQA
            19: 'http://archive.stsci.edu/pub/hlsp/phat/brick19/hlsp_phat_hst_wfc3-uvis-acs-wfc-wfc3-ir_12110-m31-b19_f275w-f336w-f475w-f814w-f110w-f160w_v2_st.fits',  # NOQA
            20: 'http://archive.stsci.edu/pub/hlsp/phat/brick20/hlsp_phat_hst_wfc3-uvis-acs-wfc-wfc3-ir_12112-m31-b20_f275w-f336w-f475w-f814w-f110w-f160w_v2_st.fits',  # NOQA
            21: 'http://archive.stsci.edu/pub/hlsp/phat/brick21/hlsp_phat_hst_wfc3-uvis-acs-wfc-wfc3-ir_12055-m31-b21_f275w-f336w-f475w-f814w-f110w-f160w_v2_st.fits',  # NOQA
            22: 'http://archive.stsci.edu/pub/hlsp/phat/brick22/hlsp_phat_hst_wfc3-uvis-acs-wfc-wfc3-ir_12076-m31-b22_f275w-f336w-f475w-f814w-f110w-f160w_v2_st.fits',  # NOQA
            23: 'http://archive.stsci.edu/pub/hlsp/phat/brick23/hlsp_phat_hst_wfc3-uvis-acs-wfc-wfc3-ir_12070-m31-b23_f275w-f336w-f475w-f814w-f110w-f160w_v2_st.fits'}  # NOQA
    url = urls[brick]
    output_path = os.path.join(os.getenv('PHATV2DIR'), os.path.basename(url))
    cmd = 'wget -c -nc -O {output} {input}'.format(output=output_path,
                                                   input=url)
    subprocess.call(cmd, shell=True)


def get_patch_info(json, brick, patch_num):
    name = '{0:02d}_{1:03d}'.format(brick, patch_num)
    for d in json:
        if d['patch'] == name:
            return d


def fit_patch(patch_info):
    """Fit a patch."""
    isoc = dict(isoc_kind='parsec_CAF09_v1.2S',
                photsys_version='yang')
    kwargs = {}
    kwargs.update(patch_info)
    kwargs['isoc_args'] = isoc
    kwargs['root_dir'] = patch_info['patch']
    pipeline = ThreeZPipeline(**kwargs)
    dataset = PatchCatalog(**patch_info)
    pipeline.fit(['oir_all'], ['oir_all'], dataset)
    pipeline.fit(['lewis'], ['lewis'], dataset)

    # TODO Reduce fits into HDF5 file object


def upload_result(result_hdf5_path, vos_dir):
    """Upload fit dataset to HDF5."""
    pass


if __name__ == '__main__':
    main()
