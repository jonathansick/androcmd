#!/usr/bin/env python
# encoding: utf-8
"""
Fit all of M31, divided into fields.

Relies on having a JSON file specifying each field. Spec is:

    {patch: <int> patch number,
     poly: <list of (ra,dec) tuples> footprint,
     brick: <int> phat brick number that photometry belongs to,
     ra0: central RA of the patch,
     dec0: central DEC of the patch,
     area: deprojected area of patch in pc^2,  # does Lewis do deproj areas?
    }

2015-05-28 - Created by Jonathan Sick
"""


def main():
    pass


if __name__ == '__main__':
    main()
