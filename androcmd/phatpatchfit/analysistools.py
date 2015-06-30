# encoding: utf-8
"""
Tools for analyzing SFR patch fits

2015-06-30 - Created by Jonathan Sick
"""

import numpy as np


def marginalize_metallicity(patch_group, fit_key):
    sfh_table = patch_group['sfh'][fit_key]
    t = np.empty(len(sfh_table), dtype=sfh_table.dtype)
    sfh_table.read_direct(t, source_sel=None, dest_sel=None)
    age_vals = np.unique(t['log(age)'])
    s = np.argsort(age_vals)
    age_vals = age_vals[s]
    A = []
    sfr = []
    for i, age_val in enumerate(age_vals):
        tt = t[t['log(age)'] == age_val]
        bin_sfr = np.sum(tt['sfr'])
        A.append(age_val)
        sfr.append(bin_sfr)
    srt = np.argsort(A)
    A = np.array(A)
    sfr = np.array(sfr)
    A = A[srt]
    sfr = sfr[srt]
    return A, sfr
