# encoding: utf-8
"""
Helper code for the baseline experiment.

2015-05-15 - Created by Jonathan Sick
"""

from starfisher.pipeline import PipelineBase
from androcmd.planes import BaselineTestPhatPlanes
from androcmd.phatpipeline import (PhatCatalog, SolarZIsocs, SolarLockfile,
                                   LewisBrickDust, PhatCrowding,
                                   ExtendedSolarIsocs, ExtendedSolarLockfile)


class SolarZPipeline(BaselineTestPhatPlanes, SolarZIsocs,
                     SolarLockfile, LewisBrickDust, PhatCrowding,
                     PipelineBase):
    """Pipeline for the baseline test with only solar metallicity."""
    def __init__(self, **kwargs):
        super(SolarZPipeline, self).__init__(**kwargs)


class ThreeZPipeline(BaselineTestPhatPlanes, ExtendedSolarIsocs,
                     ExtendedSolarLockfile, LewisBrickDust, PhatCrowding,
                     PipelineBase):
    """Pipeline for baseline test with three metallicity tracks."""
    def __init__(self, **kwargs):
        super(ThreeZPipeline, self).__init__(**kwargs)
