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
import re
from datetime import datetime

import h5py

from androcmd.phatpatchfit import PatchCatalog, ThreeZPipeline


def main():
    args = parse_args()

    # download the star catalog for this brick, if necessary
    download_brick_catalog(args.brick)

    patches = [int(v) for v in re.findall(r"[\w']+", args.patches)]

    with open(args.json_patch_path, 'r') as f:
        patch_json = json.load(f)

    for patch_num in patches:
        # fit this patch
        patch_info = get_patch_info(patch_json, args.brick, patch_num)
        result_hdf5_path = fit_patch(patch_info)

        if args.vodir is not None:
            subprocess.call('getCert', shell=True)
            upload_result(result_hdf5_path, args.vodir)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('brick', type=int,
                        help='Brick number')
    parser.add_argument('--patches',
                        help='Patch number(s) to fit in brick '
                             'can be comma-delimited')
    parser.add_argument('--json', dest='json_patch_path',
                        help='Path to patch JSON file')
    parser.add_argument('--vodir',
                        help='VOSpace directory to save results in',
                        default='vos:jonathansick/phat/patches')
    return parser.parse_args()


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
    output_path = os.path.join(os.getenv('PHATV2DATA'), os.path.basename(url))
    print "Downloading {url}".format(url)
    cmd = 'wget -c -nc -q -O {output} {input}'.format(output=output_path,
                                                      input=url)
    print "Started at", datetime.utcnow()
    if not os.path.exists(output_path):
        subprocess.call(cmd, shell=True)
    print "Ended at ", datetime.utcnow()


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

    fit_keys = ['oir_all', 'lewis']
    for fit_key in fit_keys:
        pipeline.fit(fit_key, [fit_key], dataset)

    # Output datafile
    h5path = os.path.join(os.getenv('STARFISH'), patch_info['patch'],
                          patch_info['patch'] + '.hdf5')
    hdf5 = h5py.File(h5path, mode='w')
    for k, v in patch_info.iteritems():
        hdf5.attrs[k] = v

    # Get the SFH table, making an HDF5 group
    reduce_sfh_tables(hdf5, pipeline, fit_keys)

    # Get the Hess plane of the fits
    reduce_fitted_hess_planes(hdf5, pipeline, dataset, fit_keys)

    # Save and upload hdf5 file
    hdf5.flush()

    return h5path


def reduce_sfh_tables(hdf5, pipeline, fit_keys):
    grp = hdf5.create_group('sfh')
    for fit_key in fit_keys:
        t = pipeline.fits[fit_key].solution_table(split_z=False)
        dset = grp.create_dataset(fit_key, data=t)
        dset.attrs['mean_age'] = pipeline.fits[fit_key].mean_age
    return grp


def reduce_fitted_hess_planes(hdf5, pipeline, dataset, fit_keys):
    sim_group = hdf5.create_group('sim_hess')
    obs_group = hdf5.create_group('obs_hess')
    chi_group = hdf5.create_group('chi_hess')
    diff_group = hdf5.create_group('diff_hess')

    for fit_key in fit_keys:
        sim_hess = pipeline.make_sim_hess(fit_key)
        d = _make_hess_dataset(sim_group, fit_key, sim_hess)

        obs_hess = pipeline.make_obs_hess(dataset, fit_key)
        d = _make_hess_dataset(obs_group, fit_key, obs_hess)

        diff_hess = pipeline.make_fit_diff_hess(dataset, fit_key, fit_key)
        d = _make_hess_dataset(diff_group, fit_key, diff_hess)

        chi_hess = pipeline.make_chisq_hess(dataset, fit_key, fit_key)
        chi_red = pipeline.compute_fit_chi(dataset, fit_key, fit_key,
                                           chi_hess=chi_hess)
        d = _make_hess_dataset(chi_group, fit_key, chi_hess)
        d.attrs['chi_red'] = chi_red


def _make_hess_dataset(group, fit_key, hess):
    d = group.create_dataset(fit_key, data=hess.masked_hess)
    d.attrs['origin'] = hess.origin
    d.attrs['extent'] = hess.extent
    plane = hess._plane
    d.attrs['suffix'] = plane.suffix
    d.attrs['x_mag'] = plane.x_mag
    d.attrs['y_mag'] = plane.y_mag
    d.attrs['x_span'] = plane.x_span
    d.attrs['y_span'] = plane.y_span
    d.attrs['x_label'] = plane.x_label
    d.attrs['y_label'] = plane.y_label
    d.attrs['dpix'] = plane.dpix
    return d


def upload_result(result_hdf5_path, vodir):
    """Upload fit dataset to HDF5."""
    cmd = 'vcp {0} {2}/{1}'.format(
        result_hdf5_path, os.path.basename(result_hdf5_path),
        vodir)
    print cmd
    subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    main()
