#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: staralt.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.coordinates import SkyCoord, AltAz, EarthLocation, get_body
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
    obs_location = EarthLocation.of_site('Apache Point Observatory')
    object_altaz = coords.transform_to(AltAz(obstime=date, location=obs_location))

    # Convert to float
    if isinstance(coords_string, list):
        return [alt.item() for alt in object_altaz.alt.deg]
    else:
        return object_altaz.alt.deg.item()

    # Return value for altitude
    print(f"Object's Altitude = {object_altaz.alt}")
    return object_altaz



def plot_altitude(coords_string, date_string):
    date = date_string.split()[0]
    time = Time(f"{date} 20:00:00")
    time_intervals = np.linspace(0, 12, 80) * u.hour  # 80 time steps over 12 hours
    all_times = time + time_intervals
    observer_location = EarthLocation.of_site('Apache Point Observatory')

    altaz_frame = AltAz(obstime=all_times, location=observer_location)

    plt.figure(figsize=(10, 6))

    for coord_str in coords_string:
        obj = SkyCoord(coord_str, unit=(u.hourangle, u.deg), frame='icrs')
        obj_altaz = obj.transform_to(altaz_frame)
        plt.plot(all_times.datetime, obj_altaz.alt.deg, label=f"Object: {coord_str}")

    moon = get_body('moon', time=all_times, location=observer_location)
    moon_altaz = moon.transform_to(altaz_frame)
    plt.plot(all_times.datetime, moon_altaz.alt.deg, label='Moon', color='red')

    # Final plot formatting
    plt.axhline(0, color='gray', linestyle='--', label='Horizon')
    plt.title("Altitude Plot (Object and Moon)")
    plt.xlabel("Time")
    plt.ylabel("Altitude (degrees)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.xticks(rotation=45)
    plt.savefig('altitude_plot.png')
    plt.show()

