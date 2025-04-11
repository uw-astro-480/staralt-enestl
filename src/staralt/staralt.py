#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: staralt.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy.coordinates import SkyCoord, AltAz, EarthLocation, get_body, get_sun
from astropy.coordinates import ICRS, Galactic, FK4, FK5
import astropy.units as u


def get_altitude(coords_string, date_string):
    #Takes inputs for form: (["1:33:50.89 +30:39:36.633"], "2025-4-10 23:00:00")
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


def plot_altitude():
    # Set initial variables
    coords_string = ["1:33:50.89 +30:39:36.633"]
    date = "2025-4-11"
    time = Time(F"{date} 20:00:00")
    time_intervals = np.linspace(0, 16, 80) * u.hour
    all_times = time + time_intervals
    observer_location = EarthLocation.of_site('Apache Point Observatory')

    altaz_frame = AltAz(obstime=all_times, location=observer_location)

    plt.figure(figsize=(10, 6))

    # Gets all coordinates over the time intervals
    for coord in coords_string:
        obj = SkyCoord(coord, unit=(u.hourangle, u.deg), frame='icrs')
        obj_altaz = obj.transform_to(altaz_frame)
        plt.plot(all_times.datetime, obj_altaz.alt.deg, label=f"Object: {coord}")

    # Plot moon data
    moon = get_body('moon', time=all_times, location=observer_location)
    moon_altaz = moon.transform_to(altaz_frame)
    plt.plot(all_times.datetime, moon_altaz.alt.deg, label='Moon', color='red')

    # Twilight and Sunset
    sun = get_sun(all_times)
    sun_altaz = sun.transform_to(altaz_frame)

    sunset_idx = np.where(np.diff(np.sign(sun_altaz.alt.deg)) < 0)[0]
    if len(sunset_idx) > 0:
        sunset_time = all_times[sunset_idx[0]]
        plt.axvline(sunset_time.datetime, color='yellow', linestyle='--', label='Sunset')

    twilight_idx = np.where((sun_altaz.alt.deg < -6) & (np.roll(sun_altaz.alt.deg, 1) > -6))[0]
    if len(twilight_idx) > 0:
        twilight_time = all_times[twilight_idx[0]]
        plt.axvline(twilight_time.datetime, color='purple', linestyle='--', label='Twilight')

    # LST
    #st_hours = observer_location.local_sidereal_time(all_times)
    lst_hours = all_times.sidereal_time('mean', longitude=observer_location.lon)
    lst_labels = [f"{lst:.1f}" for lst in lst_hours.hour]
    for i in range(0, len(all_times), 15):  # sparse labeling
        plt.text(all_times.datetime[i], -15, f"LST: {lst_labels[i]}h", rotation=90, fontsize=8, ha='center', color='gray')


    # Generate plot and save to folder
    plt.title("Altitude Plot (Object and Moon)")
    plt.xlabel("Time (Date-Hour)")
    plt.ylabel("Altitude (degrees)")
    plt.axhline(0, color='gray', linestyle='--', label='Horizon')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.xticks(rotation=45)
    plt.savefig('staralt.png')
    plt.show()

