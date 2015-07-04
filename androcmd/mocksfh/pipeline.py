# encoding: utf-8
"""
Pipeline to create and fit mock CMDs with StarFISH.

2015-07-03 - Created by Jonathan Sick
"""

import os

import h5py

from starfisher.testpop import TestPop


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

    def run_fit(self, fit_keys):
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
