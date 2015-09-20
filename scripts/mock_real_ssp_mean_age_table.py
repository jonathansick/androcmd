#!/usr/bin/env python
# encoding: utf-8
"""
Create a table of mean ages for the mock SSPs.

2015-09-01 - Created by Jonathan Sick
"""

# import os
import argparse
from collections import OrderedDict

import numpy as np
import h5py
from astropy.table import Table


def main():
    args = parse_args()

    hdf5 = h5py.File(args.mock_trial_set, mode='r')
    write_mean_age_table(args.output_name, hdf5['mocksfh'])


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('output_name',
                        help='Name of output table.')
    parser.add_argument('mock_trial_set',
                        help='Filename of the mock trial dataset, e.g., '
                             '`m6.fits`')
    return parser.parse_args()


def write_mean_age_table(output_path, dataset):
    sfh_list = OrderedDict([('ssp_100myr_solar', 100.),
                            ('ssp_250myr_solar', 250.),
                            ('ssp_500myr_solar', 500.),
                            ('ssp_750myr_solar', 750.),
                            ('ssp_1000myr_solar', 1000.),
                            ('ssp_1500myr_solar', 1500.),
                            ('ssp_2000myr_solar', 2000.),
                            ('ssp_3000myr_solar', 3000.)])
    plane_keys = ['lewis', 'oir_all']
    fit_mean_ages = {k: np.empty(len(sfh_list)) for k in plane_keys}
    mock_mean_ages = np.empty(len(sfh_list))
    for i, (sfh_key, mock_age_myr) in enumerate(sfh_list.items()):
        mock_mean_ages[i] \
            = dataset[sfh_key]['mock_sfh_marginal'].attrs['mean_age'][0]
        for plane_key in plane_keys:
            fit_mean_ages[plane_key][i] \
                = dataset[sfh_key]['sfh'][plane_key].attrs['mean_age'][0]
    cols = [mock_mean_ages]
    for plane_key in plane_keys:
        cols.append(fit_mean_ages[plane_key] - mock_mean_ages)
    colnames = ['Mock SSP Age (Gyr)',
                r'$\Delta A_{\texttt{ACSMS} - \mathrm{mock}}$ (Gyr)',
                r'$\Delta A_{\texttt{OIRALL} - \mathrm{mock}}$ (Gyr)']
    formats = {c: '.2f' for c in colnames}
    tbl = Table(cols, names=colnames)
    tbl.write(output_path, format='ascii.latex', formats=formats)


if __name__ == '__main__':
    main()
