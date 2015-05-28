# encoding: utf-8
"""
Tools for the pan-M31 PHAT patch fitting exercise.

2015-05-28 - Created by Jonathan Sick
"""

from collections import OrderedDict

import numpy as np

from starfisher.pipeline import PipelineBase
from starfisher.pipeline import PlaneBase
from starfisher.pipeline import ExtinctionBase
from starfisher.dust import ExtinctionDistribution

# from starfisher.pipeline import (
#     PipelineBase, IsochroneSetBase, DatasetBase, LockBase,
#     CrowdingBase, ExtinctionBase)

from androcmd.planes import make_f475w_f160w

from androcmd.phatpipeline import (PhatCrowding,
                                   ExtendedSolarIsocs, ExtendedSolarLockfile)
from androcmd.dust import mw_Av, phat_rel_extinction, LewisDustLaw


class BaselineTestPhatPlanes(PlaneBase):
    """Color plane set for the PHAT color baseline comparison test."""
    def __init__(self, **kwargs):
        self._planes = OrderedDict([
            ('oir_all', make_f475w_f160w()),
        ])
        super(BaselineTestPhatPlanes, self).__init__(**kwargs)

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


class ThreeZPipeline(BaselineTestPhatPlanes, ExtendedSolarIsocs,
                     ExtendedSolarLockfile, LewisPatchDust, PhatCrowding,
                     PipelineBase):
    """Pipeline for patch fitting with three metallicity tracks."""
    def __init__(self, **kwargs):
        super(ThreeZPipeline, self).__init__(**kwargs)
