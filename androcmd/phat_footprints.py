# encoding: utf-8
"""
Load information about PHAT footprints.

2015-06-24 - Created by Jonathan Sick
"""

import json
from pkg_resources import resource_stream, resource_exists


def load_brick_footprints():
    path = "data/phat_brick_footprints.json"
    assert resource_exists(__name__, path)
    return json.load(resource_stream(__name__, path))
