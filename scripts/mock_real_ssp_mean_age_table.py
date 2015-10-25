#!/usr/bin/env python
# encoding: utf-8
"""
Create a table of mean ages for the mock SSPs.

2015-09-01 - Created by Jonathan Sick
"""

# import os
import argparse
# from collections import OrderedDict

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
    # mock_list = [u'ssp_53myr_solar',
    #              u'ssp_100myr_solar',
    #              u'ssp_186myr_solar',
    #              u'ssp_346myr_solar',
    #              u'ssp_645myr_solar',
    #              u'ssp_1380myr_solar',
    #              u'ssp_3467myr_solar',
    #              u'ssp_5370myr_solar',
    #              u'ssp_9549myr_solar',
    #              u'tau_0.1_solar',
    #              u'tau_0.5_solar',
    #              u'tau_1.0_solar',
    #              u'tau_5.0_solar',
    #              u'tau_10.0_solar',
    #              u'tau_20.0_solar',
    #              u'tau_50.0_solar',
    #              u'tau_100.0_solar']
    mock_list = [u'ssp_53myr_solar',
                 u'ssp_100myr_solar',
                 u'ssp_186myr_solar',
                 u'ssp_346myr_solar',
                 u'ssp_645myr_solar',
                 u'ssp_1380myr_solar',
                 u'ssp_3467myr_solar',
                 u'ssp_5370myr_solar',
                 u'ssp_9549myr_solar']

    plane_keys = ['lewis', 'oir_all']
    fit_mean_ages = {k: np.empty(len(mock_list)) for k in plane_keys}
    fit_sigma_ages = {k: np.empty(len(mock_list)) for k in plane_keys}
    mock_mean_ages = np.empty(len(mock_list))
    for i, sfh_key in enumerate(mock_list):
        mock_mean_ages[i] \
            = dataset[sfh_key]['mock_sfh_marginal'].attrs['mean_age'][0]
        for plane_key in plane_keys:
            fit_mean_ages[plane_key][i] \
                = dataset[sfh_key]['sfh'][plane_key].attrs['mean_age'][0]
            fit_sigma_ages[plane_key][i] \
                = dataset[sfh_key]['sfh'][plane_key].attrs['mean_age'][1]
    cols = [mock_mean_ages]
    for plane_key in plane_keys:
        cols.append(fit_mean_ages[plane_key] - mock_mean_ages)
        cols.append(fit_sigma_ages[plane_key])
    colnames = ['Mock SSP Age (Gyr)',
                r'$\Delta A_{\texttt{ACSMS} - \mathrm{mock}}$ (Gyr)',
                r'$\sigma A_\texttt{ACSMS}$ (Gyr)',
                r'$\Delta A_{\texttt{OIRALL} - \mathrm{mock}}$ (Gyr)',
                r'$\sigma A_\texttt{OIRALL}$ (Gyr)']
    formats = {c: '.2f' for c in colnames}
    tbl = Table(cols, names=colnames)
    tbl.write(output_path, format='ascii.latex', formats=formats)


if __name__ == '__main__':
    main()
