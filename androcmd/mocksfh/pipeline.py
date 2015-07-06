# encoding: utf-8
"""
Pipeline to create and fit mock CMDs with StarFISH.

2015-07-03 - Created by Jonathan Sick
"""

import os
from collections import OrderedDict, namedtuple

import h5py

from starfisher.testpop import TestPop
from starfisher.pipeline import PipelineBase, CrowdingBase, PlaneBase
from starfisher import ColorPlane

from ..planes import make_f475w_f160w, make_lewis_ms

from ..phatpipeline import ExtendedSolarIsocs, ExtendedSolarLockfile
from ..phatpatchfit.pipeline import LewisPatchDust, AutoPhatCrowding
# from ..dust import mw_Av, phat_rel_extinction, LewisDustLaw

# from ..planes import make_f475w_f160w, make_lewis_ms

Lim = namedtuple('Lim', 'x y')


class MockFit(object):
    """Class for organizing and running a mock-SFH experiment.

    The ``sfh_factory`` should take a lockfile as the sole argument in order
    to associate output SFH amplitudes to lockfile groups.
    """
    def __init__(self, name, sfh_factory, pipeline, n_star_amp=True):
        super(MockFit, self).__init__()
        self.n_star_amp
        self.name = name
        self.sfh_factory = sfh_factory
        self.pipeline = pipeline
        self.dataset = None

    def make_dataset(self):
        sfh_amps = self.sfh_factory(self.pipeline.synth.lockfile)
        self._testpop = TestPop(self.name,
                                self.pipeline.synth,
                                sfh_amps,
                                use_lockfile=True, delta_dmod=0.,
                                n_stars=10000, fext=1., gamma=-1.35,
                                fbinary=0.5, n_star_amp=self.n_star_amp)
        self._testpop.run()
        self.dataset = self._testpop.dataset

    def run_fit(self, fit_keys, n_synth_cpu=1):
        # Ensure synth planes are prepared
        self.pipeline.run_synth(n_cpu=n_synth_cpu)

        for fit_key in fit_keys:
            self.pipeline.fit(fit_key, [fit_key], self.dataset)

        # Output datafile
        h5path = os.path.join(os.getenv('STARFISH'), self.pipeline.root_dir,
                              '{0}.hdf5'.format(self.name))
        hdf5 = h5py.File(h5path, mode='w')

        # Get the SFH table, making an HDF5 group
        self._reduce_sfh_tables(hdf5, fit_keys)

        # Get the Hess plane of the fits
        self._reduce_fitted_hess_planes(hdf5, fit_keys)

        # Save and upload hdf5 file
        hdf5.flush()

    def _reduce_sfh_tables(self, hdf5, fit_keys):
        grp = hdf5.create_group('sfh')
        for fit_key in fit_keys:
            t = self.pipeline.fits[fit_key].solution_table(split_z=False)
            dset = grp.create_dataset(fit_key, data=t)
            dset.attrs['mean_age'] = self.pipeline.fits[fit_key].mean_age
        return grp

    def _reduce_fitted_hess_planes(self, hdf5, fit_keys):
        sim_group = hdf5.create_group('sim_hess')
        obs_group = hdf5.create_group('obs_hess')
        chi_group = hdf5.create_group('chi_hess')
        diff_group = hdf5.create_group('diff_hess')

        for fit_key in fit_keys:
            sim_hess = self.pipeline.make_sim_hess(fit_key)
            d = self._make_hess_dataset(sim_group, fit_key, sim_hess)

            obs_hess = self.pipeline.make_obs_hess(self.dataset, fit_key)
            d = self._make_hess_dataset(obs_group, fit_key, obs_hess)

            diff_hess = self.pipeline.make_fit_diff_hess(self.dataset,
                                                         fit_key,
                                                         fit_key)
            d = self._make_hess_dataset(diff_group, fit_key, diff_hess)

            chi_hess = self.pipeline.make_chisq_hess(self.dataset,
                                                     fit_key, fit_key)
            chi_red = self.pipeline.compute_fit_chi(self.dataset,
                                                    fit_key, fit_key,
                                                    chi_hess=chi_hess)
            d = self._make_hess_dataset(chi_group, fit_key, chi_hess)
            d.attrs['chi_red'] = chi_red

    def _make_hess_dataset(self, group, fit_key, hess):
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


def make_f475w_f160w_28(dpix=0.05, mag_lim=36.):
    lim = Lim(x=(-0.8, 8.), y=(28., 17.5))
    plane = ColorPlane(('F475W', 'F160W'), 'F160W',
                       lim.x,
                       (min(lim.y), max(lim.y)),
                       mag_lim,
                       suffix='oira28',
                       x_label=r'$\mathrm{F475W}-\mathrm{F160W}$',
                       y_label=r'$\mathrm{F160W}$',
                       dpix=dpix)
    return plane


def make_f475w_f160w_30(dpix=0.05, mag_lim=38.):
    lim = Lim(x=(-0.8, 8.), y=(30., 17.5))
    plane = ColorPlane(('F475W', 'F160W'), 'F160W',
                       lim.x,
                       (min(lim.y), max(lim.y)),
                       mag_lim,
                       suffix='oira30',
                       x_label=r'$\mathrm{F475W}-\mathrm{F160W}$',
                       y_label=r'$\mathrm{F160W}$',
                       dpix=dpix)
    return plane


def make_f475w_f160w_32(dpix=0.05, mag_lim=40.):
    lim = Lim(x=(-0.8, 8.), y=(32., 17.5))
    plane = ColorPlane(('F475W', 'F160W'), 'F160W',
                       lim.x,
                       (min(lim.y), max(lim.y)),
                       mag_lim,
                       suffix='oira32',
                       x_label=r'$\mathrm{F475W}-\mathrm{F160W}$',
                       y_label=r'$\mathrm{F160W}$',
                       dpix=dpix)
    return plane


def make_lewis_ms_28(dpix=0.05, mag_lim=30.):
    lim = Lim(x=(-0.5, 1.), y=(28, 21.))
    plane = ColorPlane(('F475W', 'F814W'), 'F475W',
                       lim.x,
                       (min(lim.y), max(lim.y)),
                       mag_lim,
                       suffix='lws28',
                       x_label=r'$\mathrm{F475W}-\mathrm{F814W}$',
                       y_label=r'$\mathrm{F475W}$',
                       dpix=dpix)
    return plane


def make_lewis_ms_30(dpix=0.05, mag_lim=32.):
    lim = Lim(x=(-0.5, 1.), y=(30, 21.))
    plane = ColorPlane(('F475W', 'F814W'), 'F475W',
                       lim.x,
                       (min(lim.y), max(lim.y)),
                       mag_lim,
                       suffix='lws30',
                       x_label=r'$\mathrm{F475W}-\mathrm{F814W}$',
                       y_label=r'$\mathrm{F475W}$',
                       dpix=dpix)
    return plane


def make_lewis_ms_32(dpix=0.05, mag_lim=34.):
    lim = Lim(x=(-0.5, 1.), y=(32, 21.))
    plane = ColorPlane(('F475W', 'F814W'), 'F475W',
                       lim.x,
                       (min(lim.y), max(lim.y)),
                       mag_lim,
                       suffix='lws30',
                       x_label=r'$\mathrm{F475W}-\mathrm{F814W}$',
                       y_label=r'$\mathrm{F475W}$',
                       dpix=dpix)
    return plane


class DeepMockPlanes(PlaneBase):
    """Color plane set PHAT mock testing.
    """
    def __init__(self, **kwargs):
        self._planes = OrderedDict([
            ('oir_all', make_f475w_f160w()),
            ('oir_all_28', make_f475w_f160w_28()),
            ('oir_all_30', make_f475w_f160w_30()),
            ('oir_all_32', make_f475w_f160w_32()),
            ('lewis', make_lewis_ms()),
            ('lewis_28', make_lewis_ms_28()),
            ('lewis_30', make_lewis_ms_30()),
            ('lewis_32', make_lewis_ms_32()),
        ])
        super(DeepMockPlanes, self).__init__(**kwargs)

    @property
    def planes(self):
        return self._planes


class MockPlanes(PlaneBase):
    """Color plane set PHAT mock testing.

    Includes the OIR-ALL and ACS-MS planes used for the actual fitting in
    addition to OIR-ALL and ACS-MS planes that extend further down the
    luminosity function.
    """
    def __init__(self, **kwargs):
        self._planes = OrderedDict([
            ('oir_all', make_f475w_f160w()),
            ('lewis', make_lewis_ms()),
        ])
        super(MockPlanes, self).__init__(**kwargs)

    @property
    def planes(self):
        return self._planes


class RealErrorsThreeZPipeline(MockPlanes, ExtendedSolarIsocs,
                               ExtendedSolarLockfile, LewisPatchDust,
                               AutoPhatCrowding,
                               PipelineBase):
    """Pipeline for fitting with three metallicity tracks that emulates real
    fits by using PHAT crowding tables for the mock dataset.
    """
    def __init__(self, **kwargs):
        # Get patch attributes
        self.patch = kwargs.pop('patch')
        self.poly = kwargs.pop('poly')
        self.brick = kwargs.pop('brick')
        self.ra0 = kwargs.pop('ra0')
        self.dec0 = kwargs.pop('dec0')
        self.area = kwargs.pop('area')
        super(RealErrorsThreeZPipeline, self).__init__(**kwargs)


class IdealizedThreeZPipeline(DeepMockPlanes, ExtendedSolarIsocs,
                              ExtendedSolarLockfile, LewisPatchDust,
                              CrowdingBase,
                              PipelineBase):
    """Pipeline for fitting with three metallicity tracks given no crowding
    errors, and with deep CMD planes.
    """
    def __init__(self, **kwargs):
        # Get patch attributes
        self.patch = kwargs.pop('patch')
        self.poly = kwargs.pop('poly')
        self.brick = kwargs.pop('brick')
        self.ra0 = kwargs.pop('ra0')
        self.dec0 = kwargs.pop('dec0')
        self.area = kwargs.pop('area')
        super(IdealizedThreeZPipeline, self).__init__(**kwargs)
