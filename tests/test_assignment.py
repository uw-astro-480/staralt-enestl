#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: test_assignment.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)

import pathlib

from astropy.coordinates import SkyCoord
from astropy.time import Time

from staralt import get_altitude, plot_altitude


def test_get_altitude():
    """Tests the get_altitude function."""

    sirius = SkyCoord(ra=101.287155, dec=-16.716115, unit="deg")
    now = Time.now()

    altitude = get_altitude(sirius, now)

    assert altitude is not None, "Altitude should not be None"
    assert isinstance(altitude, float), "Altitude should be a float"


def test_plot_altitude():
    """Tests the plot_altitude function."""

    plot_altitude()

    output_path = pathlib.Path("./staralt.png")
    assert output_path.exists(), "Output plot file should exist"
