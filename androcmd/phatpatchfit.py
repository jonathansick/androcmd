# encoding: utf-8
"""
Tools for the pan-M31 PHAT patch fitting exercise.

2015-05-28 - Created by Jonathan Sick
"""

from collections import OrderedDict

from starfisher.pipeline import PipelineBase
from starfisher.pipeline import PlaneBase

from androcmd.planes import make_f475w_f160w

from androcmd.phatpipeline import (LewisBrickDust, PhatCrowding,
                                   ExtendedSolarIsocs, ExtendedSolarLockfile)


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


class ThreeZPipeline(BaselineTestPhatPlanes, ExtendedSolarIsocs,
                     ExtendedSolarLockfile, LewisBrickDust, PhatCrowding,
                     PipelineBase):
    """Pipeline for patch fitting with three metallicity tracks."""
    def __init__(self, **kwargs):
        super(ThreeZPipeline, self).__init__(**kwargs)
