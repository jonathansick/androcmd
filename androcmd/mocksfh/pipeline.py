# encoding: utf-8
"""
Pipeline to create and fit mock CMDs with StarFISH.

2015-07-03 - Created by Jonathan Sick
"""

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
