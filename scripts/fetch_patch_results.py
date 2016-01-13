#!/usr/bin/env python
# encoding: utf-8
"""
Get HDF5 patch results on VOSpace and combine them into a single HDF5 file.

Run as fetch_patch_results.py phat_field_starfish.hdf5  --vodir phat/fields

2015-06-16 - Created by Jonathan Sick
"""

import os
import argparse

import vos
import h5py
import numpy as np
from astropy.table import Table

from starfisher.sfh import estimate_mean_age, marginalize_sfh_metallicity

from androcmd.phatpatchfit import compute_patch_gal_coords


def main():
    args = parse_args()

    dataset = h5py.File(args.output_hdf5, 'a')
    dataset_patches = dataset.require_group('patches')

    vosdir = os.path.join('vos:jonathansick', args.vodir)
    c = vos.Client('vos:jonathansick')
    file_names = c.listdir(vosdir)
    for name in file_names:
        print name
        if not name.endswith('.hdf5'):
            continue

        patch_name = os.path.splitext(os.path.basename(name))[0]
        if patch_name in dataset_patches:
            print "{0} already in dataset".format(patch_name)
            continue

        try:
            c.copy(os.path.join(vosdir, name), name)
        except:
            print "Could not download {0}".format(name)
            continue

        # open the patch's HDF5 file
        try:
            patch_data = h5py.File(name, 'r')
        except IOError:
            continue
        patch_serial = patch_data.attrs['patch']
        patch_data.copy(patch_data,
                        dataset_patches,
                        name=patch_serial)
        patch_data.close()
        os.remove(name)

    validate_coords(dataset_patches)

    reduce_sfh_table(dataset, dataset_patches)

    dataset.flush()
    dataset.close()


def parse_args():
    parser = argparse.ArgumentParser(
        description="fetch_patch_results.py phat_field_starfish.hdf5 "
                    "--vodir phat/fields")
    parser.add_argument('output_hdf5')
    parser.add_argument('--vodir', default='phat/patches')
    return parser.parse_args()


def validate_coords(patch_data):
    for patch_name, group in patch_data.items():
        print patch_name, group
        if 'r_kpc' not in group.attrs.keys():
            print "Adding galaxy coords for {0}".format(patch_name)
            r_kpc, phi = compute_patch_gal_coords(group.attrs['ra0'],
                                                  group.attrs['dec0'])
            group.attrs['r_kpc'] = r_kpc
            group.attrs['phi'] = phi


def reduce_sfh_table(dataset, patches, fit_keys=None):
    patch_names = []
    ra = []
    dec = []
    r_kpc = []
    phi = []
    mean_ages = []
    mean_age_errs = []
    ages_25 = []
    ages_75 = []
    chisqs = []
    for patch_name, patch_group in patches.items():
        patch_names.append(patch_name)
        r_kpc.append(patch_group.attrs['r_kpc'])
        phi.append(patch_group.attrs['phi'])
        ra.append(patch_group.attrs['ra0'])
        dec.append(patch_group.attrs['dec0'])

        if fit_keys is None:
            fit_keys = patch_group['sfh'].keys()

        # redo the mean age estimation for each SFH; the SFH may have been
        # persisted with an incorrect mean age estimation algo
        for fit_key in fit_keys:
            sfh_table = patch_group['sfh'][fit_key]
            bin_age = sfh_table['log(age)'] ** 10. / 1e9
            mass = sfh_table['mass']
            mass_positive_sigma = sfh_table['mass_pos_err']
            mass_negative_sigma = sfh_table['mass_neg_err']
            mean_age = estimate_mean_age(
                bin_age, mass,
                mass_positive_sigma=mass_positive_sigma,
                mass_negative_sigma=mass_negative_sigma)
            sfh_table.attrs['mean_age'] = mean_age

        mean_ages.append([patch_group['sfh'][fit_key].attrs['mean_age'][0]
                          for fit_key in fit_keys])
        mean_age_errs.append([patch_group['sfh'][fit_key].attrs['mean_age'][1]
                              for fit_key in fit_keys])
        chisqs.append([patch_group['chi_hess'][fit_key].attrs['chi_red']
                       for fit_key in fit_keys])

        # compute a marginalized SFH and persist it to the dataset
        marginal_sfh_group = patch_group.create_group('sfh_marginal')
        for fit_key in fit_keys:
            sfh_table = Table(np.array(patch_group['sfh'][fit_key]))
            marginalized_sfh_table = marginalize_sfh_metallicity(sfh_table)
            marginal_sfh_group.create_dataset(
                fit_key,
                data=np.array(marginalized_sfh_table))

        # Compute the 25th and 75th percentils of cumulative SF.
        _25 = []
        _75 = []
        for fit_key in fit_keys:
            t = patch_group['sfh_marginal'][fit_key]
            age_gyr = t['log(age)'] ** 10. / 1e9
            srt = np.argsort(age_gyr)
            age_gyr = age_gyr[srt]
            mass = t['mass'][srt]
            fractional_mass = np.cumsum(mass) / mass.sum() * 100.
            result = np.interp([25., 75.],
                               fractional_mass,
                               age_gyr)
            q25 = result[0]
            q75 = result[1]
            _25.append(q25)
            _75.append(q75)
        ages_25.append(_25)
        ages_75.append(_75)

    # Build a record array
    age_fmt = 'mean_age_{0}'
    age_err_fmt = 'mean_age_err_{0}'
    age_25_fmt = 'age_25_{0}'
    age_75_fmt = 'age_75_{0}'
    chi_fmt = 'chi_red_{0}'
    dtype = [('name', 'S40'), ('r_kpc', float), ('phi', float)] \
        + [('ra', float), ('dec', float)] \
        + [(age_fmt.format(n), float) for n in fit_keys] \
        + [(age_err_fmt.format(n), float) for n in fit_keys] \
        + [(age_25_fmt.format(n), float) for n in fit_keys] \
        + [(age_75_fmt.format(n), float) for n in fit_keys] \
        + [(chi_fmt.format(n), float) for n in fit_keys]
    n = len(patch_names)
    sfh_table = np.empty(n, dtype=np.dtype(dtype))
    sfh_table['name'][:] = patch_names
    sfh_table['r_kpc'][:] = r_kpc
    sfh_table['phi'][:] = phi
    sfh_table['ra'][:] = ra
    sfh_table['dec'][:] = dec
    for i, fit_key in enumerate(fit_keys):
        sfh_table[age_fmt.format(fit_key)][:] = [v[i] for v in mean_ages]
        sfh_table[age_err_fmt.format(fit_key)][:] = [v[i]
                                                     for v in mean_age_errs]
        sfh_table[age_25_fmt.format(fit_key)][:] = [v[i] for v in ages_25]
        sfh_table[age_75_fmt.format(fit_key)][:] = [v[i] for v in ages_75]
        sfh_table[age_75_fmt.format(fit_key)][:] = [v[i] for v in ages_75]
        sfh_table[chi_fmt.format(fit_key)][:] = [v[i] for v in chisqs]

    if 'sfh_table' in dataset.keys():
        del dataset['sfh_table']
    dataset.create_dataset('sfh_table', data=sfh_table)


if __name__ == '__main__':
    main()
