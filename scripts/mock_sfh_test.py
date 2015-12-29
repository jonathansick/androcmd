#!/usr/bin/env python
# encoding: utf-8
"""
Run mock star formation history fitting tests.

2015-07-04 - Created by Jonathan Sick
"""

import os
import json
import argparse
from functools import partial
from multiprocessing import Pool

import numpy as np
import h5py

from androcmd.mocksfh.pipeline import (RealErrorsThreeZPipeline,
                                       IdealizedThreeZPipeline)
from androcmd.mocksfh.pipeline import MockFit


def main():
    args = parse_args()

    with open(args.json_patch_path, 'r') as f:
        patch_json = json.load(f)

    patch_info = get_patch_info(patch_json, args.brick, args.patch)

    # setup the pipeline, including running synth
    isoc = dict(isoc_kind='parsec_CAF09_v1.2S',
                photsys_version='yang')
    kwargs = {}
    kwargs.update(patch_info)
    kwargs['isoc_args'] = isoc
    kwargs['root_dir'] = args.name
    kwargs['n_synth_cpu'] = args.n_synth_cpu
    P = PIPELINES[args.pipeline]
    p = P(**kwargs)

    h5path = os.path.join(os.getenv('STARFISH'), p.root_dir,
                          '{0}.hdf5'.format(args.name))
    hdf5 = h5py.File(h5path, mode='a')
    exp_group = hdf5.require_group('mocksfh')

    mocks = {}
    fit_args = []
    for i, sfh_name in enumerate(args.sfh_names):
        factory = SFH_FACTORIES[sfh_name]
        mockfit = MockFit(args.name, factory, p, n_star_amp=True)
        mockfit.make_dataset()
        mocks[sfh_name] = mockfit
        if sfh_name not in exp_group.keys():
            # do not repeat computations
            fit_args.append((sfh_name, mockfit, args.fit, i))

    print "Fitting {0:d} star formation histories".format(len(fit_args))

    if args.n_synth_cpu > 1:
        pool = Pool(processes=args.n_synth_cpu)
        M = pool.map
    else:
        M = map
    result = M(_run_fit, fit_args)
    for product in result:
        sfh_name, mockfit = product
        print sfh_name
        if sfh_name in exp_group.keys():
            del exp_group[sfh_name]
            hdf5.flush()
        sfh_group = exp_group.create_group(sfh_name)
        mockfit.persist_fit_to_hdf5(sfh_group)
        print "Persisted {0} to HDF5".format(sfh_name)

    hdf5.flush()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('name',
                        help='Test name')
    parser.add_argument('--sfh',
                        dest='sfh_names',
                        type=str, nargs='*',
                        choices=SFH_FACTORIES.keys(),
                        default=SFH_FACTORIES.keys())
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
                        default=['lewis', 'oir_all'],
                        help='Names of fitted planes')
    return parser.parse_args()


def ssp_solar(lockfile, age_myr=100):
    """Build a SFR table for this lockfile given an SSP of the given age."""
    ages = 10 ** (lockfile.group_logages - 6.)
    Zs = lockfile.group_metallicities
    n_groups = len(ages)

    # select the single group that matches the description
    i = np.argmin(np.hypot(ages - age_myr, Zs - 0.019))

    sfhs = np.zeros(n_groups, dtype=np.float)
    sfhs[i] = 1.
    return sfhs


SFH_FACTORIES = dict()

# These ages match the ExtendedSolarIsocs isochrone set
# so that SSP ages fall correctly into isochrone bins
for myr in (53.7, 100.0, 186.2, 346.7, 645.7, 1380.4, 3467.4, 5370.3, 9549.9):
    key = 'ssp_{0:d}myr_solar'.format(int(myr))
    SFH_FACTORIES[key] = partial(ssp_solar, age_myr=myr)


def tau_solar(lockfile, tau=1., tform=10.):
    """Build a SFR table for this lockfile given a declining exponential SFH.
    """
    ages_gyr = 10 ** (lockfile.group_logages - 9.)
    Zs = lockfile.group_metallicities
    dt = lockfile.group_dt / 1e9  # time spans for each group, in Gyr
    n_groups = len(ages_gyr)

    sfhs = np.zeros(n_groups, dtype=np.float)
    z_values = np.unique(Zs)
    for Z in z_values:
        s = np.where(Zs == Z)[0]
        A = 13.7 - ages_gyr[s]  # Gyr since big bang
        age_tform = 13.7 - tform  # Gyr since big bang
        sfr = np.exp(-(A - age_tform) / tau)
        # Normalize star formation
        sfr[A < age_tform] = 0.
        sfhs[s] = sfr

    # Normalize
    sfhs = sfhs / np.sum(sfhs * dt)

    return sfhs


for tau in (0.1, 0.5, 1., 5., 10., 20., 50., 100.):
    key = 'tau_{0:.1f}_solar'.format(tau)
    SFH_FACTORIES[key] = partial(tau_solar, tau=tau, tform=11.)


PIPELINES = {
    'realistic': RealErrorsThreeZPipeline,
    'ideal': IdealizedThreeZPipeline,
}


def get_patch_info(json, brick, patch_num):
    name = '{0:02d}_{1:03d}'.format(brick, patch_num)
    for d in json:
        if d['patch'] == name:
            return d


def _run_fit(args):
    sfh_name, mockfit, fit_keys, index = args
    print '_run_fit', sfh_name, fit_keys, index
    mockfit.run_fit(fit_keys, index)
    return sfh_name, mockfit


if __name__ == '__main__':
    main()
