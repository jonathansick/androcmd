# encoding: utf-8
"""
Code for fitting and plotting star formation histories in PHAT fields

2015-06-30 - Created by Jonathan Sick
"""

from .pipeline import (PatchFitPlanes, LewisPatchDust, AutoPhatCrowding,
                       ThreeZPipeline, PatchCatalog, build_patches,
                       build_field_patches, compute_patch_gal_coords,
                       galaxy_coords, load_brick_footprints,
                       load_field_patches, load_field_footprints)
from .analysistools import marginalize_metallicity
from .galexmap import load_galex_map, setup_galex_axes, plot_patch_footprints
from .planecompmaps import setup_plane_comp_axes
from .sfrplots import get_scaled_sfr_values, scale_sfr, lin_scale_sfr, \
    SFR_LABEL, LIN_SFR_LABEL
