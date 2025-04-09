#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: staralt.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


def get_altitude(coords, date):
    """Returns the altitude of the object at the given date.

    Parameters
    ----------
    coords
        An astropy SkyCoord object with the coordinates of the object.
    date
        An astropy Time object with the date-time of the observation.

    Returns
    -------
    altitude
        The altitude of the object in degrees, as a float.

    """

    raise NotImplementedError("This function is not implemented yet.")


def plot_altitude():
    """Plots the altitude of the object over time, for a certain night."""

    raise NotImplementedError("This function is not implemented yet.")
