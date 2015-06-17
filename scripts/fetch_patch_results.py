#!/usr/bin/env python
# encoding: utf-8
"""
Get HDF5 patch results on VOSpace and combine them into a single HDF5 file.

2015-06-16 - Created by Jonathan Sick
"""

import os
import argparse

import vos
import h5py

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

    dataset.flush()
    dataset.close()


def parse_args():
    parser = argparse.ArgumentParser()
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


if __name__ == '__main__':
    main()
