# encoding: utf-8
"""
Tools for the pan-M31 PHAT patch fitting exercise.

2015-05-28 - Created by Jonathan Sick
"""

import os
from collections import OrderedDict

from matplotlib.path import Path

import numpy as np
from astropy.table import Table
from astropy.io import fits
import astropy
from astropy.coordinates import Distance, Angle, SkyCoord
import astropy.units as u


from m31hst import phat_v2_phot_path, phat_brick_path
from m31hst.phatast import PhatAstTable

from starfisher.pipeline import PipelineBase
from starfisher.pipeline import PlaneBase
from starfisher.pipeline import ExtinctionBase
from starfisher.pipeline import DatasetBase
from starfisher.pipeline import CrowdingBase
from starfisher.dust import ExtinctionDistribution
from starfisher import ExtantCrowdingTable

# from starfisher.pipeline import (
#     PipelineBase, IsochroneSetBase, DatasetBase, LockBase,
#     CrowdingBase, ExtinctionBase)

from androcmd.planes import make_f475w_f160w, make_lewis_ms

from androcmd.phatpipeline import ExtendedSolarIsocs, ExtendedSolarLockfile
from androcmd.dust import mw_Av, phat_rel_extinction, LewisDustLaw


class PatchFitPlanes(PlaneBase):
    """Color plane set for the PHAT color baseline comparison test."""
    def __init__(self, **kwargs):
        self._planes = OrderedDict([
            ('oir_all', make_f475w_f160w()),
            ('lewis', make_lewis_ms()),
        ])
        super(PatchFitPlanes, self).__init__(**kwargs)

    @property
    def planes(self):
        return self._planes


class LewisPatchDust(ExtinctionBase):
    """Mixin for a uniform dust distribution fitted from Lewis et al 15 Fig 17.

    The maximum extinction is estimated from a Draine et al 2015 dust map.
    Requires that the brick be known.
    """
    def __init__(self, **kwargs):
        # self.brick = kwargs.pop('brick', 23)
        super(LewisPatchDust, self).__init__(**kwargs)

    def build_extinction(self):
        """Young and old dust at equal here."""
        lewis = LewisDustLaw()
        max_av = lewis.estimate_mean_extinction(self.poly)
        av = np.random.uniform(low=mw_Av(),
                               high=max_av,
                               size=1000)

        self.young_av = ExtinctionDistribution()
        self.young_av.set_samples(av)

        self.old_av = ExtinctionDistribution()
        self.old_av.set_samples(av)

        self.rel_extinction = phat_rel_extinction()


class AutoPhatCrowding(CrowdingBase):
    """Use crowding from the PHAT AST fields, automatically selecting
    the AST field that best matches the location of the patch being
    considered.
    """
    def __init__(self, **kwargs):
        super(AutoPhatCrowding, self).__init__(**kwargs)

    def build_crowding(self):
        crowd_path = os.path.join(self.synth_dir, "crowding.dat")
        full_crowd_path = os.path.join(os.getenv('STARFISH'), crowd_path)

        patch_coord = SkyCoord(ra=self.ra0 * u.degree,
                               dec=self.dec0 * u.degree)

        tbl = PhatAstTable()
        centers = np.array([f['center'] for f in tbl.fields])

        ast_coords = SkyCoord(ra=centers[:, 0] * u.degree,
                              dec=centers[:, 1] * u.degree,
                              frame='icrs')
        dists = ast_coords.separation(patch_coord)
        best_ast_field = np.argmin(dists)
        # FIXME actually want to compute different in galactocentric radii?

        tbl.write_crowdfile_for_field(full_crowd_path,
                                      best_ast_field,
                                      bands=self.bands)
        self.crowd = ExtantCrowdingTable(crowd_path)

    def mask_planes(self):
        """Mask each CMD plane based on the incomplete or empty regions of
        the PHAT artificial star testing projected into the Hess plane.

        This hook is called automatically by the base pipeline before
        synth is run.
        """
        # FIXME note that AST field 0 *is always* used
        print "Using PhatCrowding.mask_planes"
        ast = PhatAstTable()
        for key, plane in self.planes.iteritems():
            band = plane.y_mag  # FIXME assumes CMD; only 1 y axis mag.
            hess, x_grid, y_grid = ast.completeness_hess(
                0, band,
                plane.x_mag, plane.y_mag,
                plane.xlim, plane.ylim, 0.5)
            yidx, xidx = np.where(hess < 0.5)  # mask less than 50% complete
            for yi, xi in zip(yidx, xidx):
                plane.mask_region((x_grid[xi], x_grid[xi + 1]),
                                  (y_grid[yi], y_grid[yi + 1]))
            yidx, xidx = np.where(~np.isfinite(hess))  # mask empty AST
            for yi, xi in zip(yidx, xidx):
                plane.mask_region((x_grid[xi], x_grid[xi + 1]),
                                  (y_grid[yi], y_grid[yi + 1]))


class ThreeZPipeline(PatchFitPlanes, ExtendedSolarIsocs,
                     ExtendedSolarLockfile, LewisPatchDust, AutoPhatCrowding,
                     PipelineBase):
    """Pipeline for patch fitting with three metallicity tracks."""
    def __init__(self, **kwargs):
        # Get patch attributes
        self.patch = kwargs.pop('patch')
        self.poly = kwargs.pop('poly')
        self.brick = kwargs.pop('brick')
        self.ra0 = kwargs.pop('ra0')
        self.dec0 = kwargs.pop('dec0')
        self.area = kwargs.pop('area')
        super(ThreeZPipeline, self).__init__(**kwargs)


class PatchCatalog(DatasetBase):
    """Mixin for PHAT photometry data.

    Photometry is lazy loaded to it is efficient to rebuild the pipeline
    object.
    """
    def __init__(self, **kwargs):
        self._phat_data = None
        self.patch = kwargs.pop('patch')
        self.poly = kwargs.pop('poly')
        self.brick = kwargs.pop('brick')
        self.ra0 = kwargs.pop('ra0')
        self.dec0 = kwargs.pop('dec0')
        self.area = kwargs.pop('area')
        self.area_proj = kwargs.pop('area_proj')
        self.r_kpc = kwargs.pop('r_kpc')
        self.phi = kwargs.pop('phi')
        print "uncaught PatchCatalog keywords", kwargs
        super(PatchCatalog, self).__init__(**kwargs)

    def _load_phat_data(self):
        self._phat_data = Table.read(phat_v2_phot_path(self.brick),
                                     format='fits')
        phat_bands = ('F475W', 'F814W', 'F275W', 'F336W', 'F110W', 'F160W')
        # Normalize bandpass names
        for band in phat_bands:
            old_name = "_".join((band.lower(), 'vega'))
            self._phat_data.rename_column(old_name, band)

        # Spatially select stars in the patch footprint
        self._phat_data = self._spatial_selection(self._phat_data)

    def _spatial_selection(self, dataset):
        data_ra = dataset['ra']
        data_dec = dataset['dec']
        xy = np.vstack((data_ra, data_dec)).T
        print "xy.shape", xy.shape

        polygon = Path(self.poly, closed=False)
        inside = polygon.contains_points(xy)
        return dataset[inside == True]  # noqa

    def get_phot(self, band):
        if self._phat_data is None:
            self._load_phat_data()  # lazy loading

        if not isinstance(band, basestring):
            band1, band2 = band
            return self._phat_data[band1] - self._phat_data[band2]
        else:
            return self._phat_data[band]

    def _select_gst(self, x_band, y_band):
        mags = []
        if not isinstance(x_band, basestring):
            mags.extend(x_band)
        else:
            mags.append(x_band)

        if not isinstance(y_band, basestring):
            mags.extend(y_band)
        else:
            mags.append(y_band)
        mags = list(set(mags))

        gsts = []
        for band in mags:
            key = '{0}_gst'.format(band.lower())
            gsts.append(self._phat_data[key])
        gsts_array = np.vstack(gsts).T
        gsts = np.all(gsts_array, axis=1)
        return np.where(gsts == True)[0]  # NOQA

    def write_phot(self, x_band, y_band, data_root, suffix):
        """Only good (GST=1 in all relevant bands) photometry is written."""
        if self._phat_data is None:
            self._load_phat_data()  # lazy loading

        x = self.get_phot(x_band)
        y = self.get_phot(y_band)
        gst_sel = self._select_gst(x_band, y_band)

        phot_dtype = np.dtype([('x', np.float), ('y', np.float)])
        photdata = np.empty(len(gst_sel), dtype=phot_dtype)
        photdata['x'][:] = x[gst_sel]
        photdata['y'][:] = y[gst_sel]

        path = data_root + suffix
        full_path = os.path.join(os.getenv('STARFISH'), path)
        fit_dir = os.path.dirname(full_path)
        if not os.path.exists(fit_dir):
            os.makedirs(fit_dir)
        np.savetxt(full_path, photdata, delimiter=' ', fmt='%.4f')

    @property
    def polygon(self):
        """Polygon bounding box around the dataset."""
        return self.poly


def build_patches(brick, proj_size=100):
    """Build patches from a brick.

    brick : int
        PHAT brick number.
    proj_size : float
        Size of a patch, in projected side length, in parsecs.
    """
    patches = []

    # Get the patch FITS file
    path = phat_brick_path(brick, 'F814W')
    with fits.open(path) as f:
        header = fits.getheader(f, 0)
        wcs = astropy.wcs.WCS(header)

    # degree per pixel
    pixel_scale = np.mean(
        np.sqrt(astropy.wcs.utils.proj_plane_pixel_scales(wcs) ** 2.))

    theta = np.rad2deg(np.arctan(proj_size / (785. * 10. ** 3.)))

    n_pix_side_min = theta / pixel_scale

    # Number of boxes, in each dimension so each box is *at least*
    # proj_size on each size.
    nx = np.floor(header['NAXIS1'] / n_pix_side_min)
    ny = np.floor(header['NAXIS2'] / n_pix_side_min)

    x_edges = np.linspace(0, header['NAXIS1'], num=nx + 1,
                          endpoint=True, dtype=int)
    y_edges = np.linspace(0, header['NAXIS2'], num=ny + 1,
                          endpoint=True, dtype=int)

    patch_num = 1
    for i in xrange(nx):
        for j in xrange(ny):
            x1 = x_edges[i]
            x2 = x_edges[i + 1]
            y1 = y_edges[i]
            y2 = y_edges[i + 1]
            xy_verts = np.array([[x1, y1],
                                 [x1, y2],
                                 [x2, y2],
                                 [x2, y1]])
            radec = wcs.all_pix2world(xy_verts, 0)
            # projected arcsec^2
            area_proj = (y2 - y1) * (x2 - x1) * (pixel_scale * 3600.) ** 2.
            # area in pc^2, de-projected
            pc_per_arcsec = 785. * 10. ** 3. * np.tan(np.deg2rad(1. / 3600))
            area = area_proj * (pc_per_arcsec) ** 2. / np.cos(77 * np.pi / 180.)  # NOQA
            ra0 = radec[:, 0].mean()
            dec0 = radec[:, 1].mean()
            r_kpc, phi = compute_patch_gal_coords(ra0, dec0)
            patch = {'patch': '{0:02d}_{1:03d}'.format(brick, patch_num),
                     'brick': brick,
                     'poly': radec.tolist(),
                     'ra0': ra0,
                     'dec0': dec0,
                     'r_kpc': r_kpc,
                     'phi': phi,
                     'area_proj': area_proj,
                     'area': area}
            patches.append(patch)
            patch_num += 1

    return patches


def compute_patch_gal_coords(ra, dec):
    """Compute the galactocentric coordinates (radius in kpc; PA in degrees)
    for a patch given its central coordinate.
    """
    d0 = Distance(785, unit=u.kpc)
    pa0 = Angle('37d42m54s')
    incl0 = Angle('77.5d')
    coord0 = SkyCoord('00h42m44.33s +41d16m07.5s')
    coord = SkyCoord(ra * u.degree, dec * u.degree)
    r, phi, r_sky = galaxy_coords(coord, coord0, pa0, incl0, d0)
    return r.kpc, phi.degree


def galaxy_coords(coord, glx_ctr, glx_PA, glx_incl, glx_dist):
    """Computes deprojected galactocentric distance.
    Inspired by: http://idl-moustakas.googlecode.com/svn-history/
        r560/trunk/impro/hiiregions/im_hiiregion_deproject.pro
    Parameters
    ----------
    coord : :class:`astropy.coordinates.SkyCoord`
        Coordinate of points to compute galactocentric distance for.
        Can be either a single coordinate, or array of coordinates.
    glx_ctr : :class:`astropy.coordinates.SkyCoord`
        Galaxy center.
    glx_PA : :class:`astropy.coordinates.Angle`
        Position angle of galaxy disk.
    glx_incl : :class:`astropy.coordinates.Angle`
        Inclination angle of the galaxy disk.
    glx_dist : :class:`astropy.coordinates.Distance`
        Distance to galaxy.
    Returns
    -------
    obj_dist : class:`astropy.coordinates.Distance`
        Galactocentric distance(s) for coordinate point(s).
    """
    # distance from coord to glx centre
    sky_radius = glx_ctr.separation(coord)
    avg_dec = 0.5 * (glx_ctr.dec + coord.dec).radian
    x = (glx_ctr.ra - coord.ra) * np.cos(avg_dec)
    y = glx_ctr.dec - coord.dec
    # azimuthal angle from coord to glx  -- not completely happy with this
    phi = glx_PA - Angle('90d') \
        + Angle(np.arctan2(y.arcsec, x.arcsec), unit=u.rad)

    # convert to coordinates in rotated frame, where y-axis is galaxy major
    # ax; have to convert to arcmin b/c can't do sqrt(x^2+y^2) when x and y
    # are angles
    xp = (sky_radius * np.cos(phi.radian)).arcmin
    yp = (sky_radius * np.sin(phi.radian)).arcmin

    # de-project
    ypp = yp / np.cos(glx_incl.radian)
    obj_radius = np.sqrt(xp ** 2 + ypp ** 2)  # in arcmin
    obj_dist = Distance(Angle(obj_radius, unit=u.arcmin).radian * glx_dist,
                        unit=glx_dist.unit)

    # Computing PA in disk
    # negative sign needed to get correct orientation from major axis
    obj_phi = Angle(np.arctan2(ypp, -xp), unit=u.rad)
    # FIXME only for scalars
    if obj_phi < 0.:
        obj_phi = Angle(2. * np.pi, unit=u.rad) + obj_phi

    return obj_dist, obj_phi, sky_radius.arcsec
