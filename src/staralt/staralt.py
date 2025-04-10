#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: staralt.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord, AltAz, EarthLocation
from astropy.coordinates import ICRS, Galactic, FK4, FK5
import astropy.units as u


def get_altitude(coords_string, date_string):
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
    # Takes inputs for form: 
    # (["1:33:50.89 +30:39:36.633", "5:14:32.3 -8:12:5.1"], "2025-4-10 23:00:00")
    
    #Transforming coords_string and date_string into SkyCoords and astropy Time variables
    coords = SkyCoord(coords_string, unit=(u.hourangle, u.deg), frame='icrs')
    date = Time(date_string)

    # Get altitude and azimuth
    obs_location = EarthLocation(lat=41.3*u.deg, lon=-74*u.deg, height=0*u.m)
    object_altaz = coords.transform_to(AltAz(obstime=date, location=obs_location))

    # Convert to float
    if isinstance(coords_string, list):
        return [alt.item() for alt in object_altaz.alt.deg]
    else:
        return object_altaz.alt.deg.item()

    # Return value for altitude
    print(f"Object's Altitude = {object_altaz.alt}")
    return object_altaz



def plot_altitude():
    """Plots the altitude of the object over time, for a certain night."""

    raise NotImplementedError("This function is not implemented yet.")
