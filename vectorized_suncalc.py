import numpy as np
import datetime
import math
from suncalc_old import julianday, refraction_at_zenith, jday_to_jcentury, geom_mean_anomaly_sun
from suncalc_old import sun_declination, SunDirection, SUN_APPARENT_RADIUS, geom_mean_long_sun
import pytz

'''
The code in this repository is adopted from code in the suncalcPy library. 
All functions were vectorized and can be performed millions of times faster. 
'''

def get_sunrise_vec(latitude, longitude, d0):
    return time_of_transit_vec(
        latitude_vec=latitude,
        longitude_vec=longitude,
        d0=d0,
        direction=SunDirection.RISING,
        zenith=90.0 + SUN_APPARENT_RADIUS
    )

def get_sunset_vec(latitude, longitude, d0):
    return time_of_transit_vec(
        latitude_vec=latitude,
        longitude_vec=longitude,
        d0=d0,
        direction=SunDirection.SETTING,
        zenith=90.0 + SUN_APPARENT_RADIUS
    )

def time_of_transit_vec(
        latitude_vec,
        longitude_vec,
        d0: datetime.date,
        zenith: float,
        direction: SunDirection,
) -> datetime.datetime:
    adjustment_for_elevation = 0
    adjustment_for_refraction = refraction_at_zenith(zenith + adjustment_for_elevation)

    jd = julianday(d0)
    t = jday_to_jcentury(jd)
    solarDec = sun_declination(t)

    hourangle = hour_angle_vec(
        latitude_vec,
        solarDec,
        zenith + adjustment_for_elevation - adjustment_for_refraction,
        direction,
    )

    delta = -longitude_vec - 180 * hourangle / math.pi
    timeDiff = 4.0 * delta
    timeUTC = 720.0 + timeDiff - eq_of_time_vec(t)

    t = jday_to_jcentury_vec(jcentury_to_jday_vec(t) + timeUTC / 1440.0)
    solarDec = sun_declination_vec(t)

    hourangle = hour_angle_vec(
        latitude_vec,
        solarDec,
        zenith + adjustment_for_elevation + adjustment_for_refraction,
        direction,
    )

    delta = -longitude_vec - 180 * hourangle / math.pi
    timeDiff = 4.0 * delta
    timeUTC = 720 + timeDiff - eq_of_time_vec(t)

    td = np.array(minutes_to_timedelta_vec(timeUTC))
    dt = datetime.datetime(d0.year, d0.month, d0.day) + td

    dt = [pytz.utc.localize(d) for d in dt]  # pylint: disable=E1120
    return dt


def jcentury_to_jday_vec(juliancentury_vec) -> np.array:
    """Convert a Julian Century number to a Julian Day"""
    return (juliancentury_vec * 36525.0) + 2451545.0


def jday_to_jcentury_vec(julianday_vec) -> float:
    """Convert a Julian Day number to a Julian Century"""
    return (julianday_vec - 2451545.0) / 36525.0


def minutes_to_timedelta_vec(minutes) -> datetime.timedelta:
    """Convert a floating point number of minutes to a :class:`~datetime.timedelta`"""
    d = (minutes / 1440).astype(int)
    minutes = minutes - (d * 1440)
    seconds = minutes * 60
    s = (seconds).astype(int)
    sfrac = seconds - s
    us = (sfrac * 1_000_000).astype(int)
    arr = np.column_stack([d, s, us])

    tds = []
    for i, a in enumerate(arr.tolist()):
        try:
            tds.append(datetime.timedelta(days=a[0], seconds=a[1], microseconds=a[2]))
        except:
            tds.append(datetime.timedelta(days=-99))#np.nan)
    return tds


def mean_obliquity_of_ecliptic_vec(juliancentury_vec) -> float:
    seconds = 21.448 - juliancentury_vec * (
            46.815 + juliancentury_vec * (0.00059 - juliancentury_vec * (0.001813))
    )
    return 23.0 + (26.0 + (seconds / 60.0)) / 60.0


def sun_declination_vec(juliancentury_vec) -> float:
    """Calculate the sun's declination"""
    e = obliquity_correction_vec(juliancentury_vec)
    lambd = sun_apparent_long_vec(juliancentury_vec)

    sint = np.sin(math.pi * e / 180) * np.sin(math.pi * lambd / 180)
    return np.arcsin(sint) * 180 / math.pi


def obliquity_correction_vec(juliancentury_vec) -> float:
    e0 = mean_obliquity_of_ecliptic_vec(juliancentury_vec)

    omega = 125.04 - 1934.136 * juliancentury_vec
    return e0 + 0.00256 * np.cos(math.pi * omega / 180)


def sun_apparent_long_vec(juliancentury: float) -> float:
    true_long = sun_true_long_vec(juliancentury)

    omega = 125.04 - 1934.136 * juliancentury
    return true_long - 0.00569 - 0.00478 * np.sin(math.pi / 180 * omega)


def sun_true_long_vec(juliancentury: float) -> float:
    """Calculate the sun's true longitude"""
    l0 = geom_mean_long_sun(juliancentury)
    c = sun_eq_of_center_vec(juliancentury)

    return l0 + c


def sun_eq_of_center_vec(juliancentury: float) -> float:
    """Calculate the equation of the center of the sun"""
    m = geom_mean_anomaly_sun(juliancentury)

    mrad = math.pi / 180 * m
    sinm = np.sin(mrad)
    sin2m = np.sin(mrad + mrad)
    sin3m = np.sin(mrad + mrad + mrad)

    c = (
            sinm * (1.914602 - juliancentury * (0.004817 + 0.000014 * juliancentury))
            + sin2m * (0.019993 - 0.000101 * juliancentury)
            + sin3m * 0.000289
    )
    return c

def eq_of_time_vec(juliancentury: float) -> float:
    l0 = geom_mean_long_sun(juliancentury)
    e = eccentric_location_earth_orbit_vec(juliancentury)
    m = geom_mean_anomaly_sun(juliancentury)

    y = var_y_vec(juliancentury)

    sin2l0 = np.sin(2.0 * math.pi / 180 * (l0))
    sinm = np.sin(math.pi / 180 * m)
    cos2l0 = np.cos(2.0 * math.pi / 180 * (l0))
    sin4l0 = np.sin(4.0 * math.pi / 180 * (l0))
    sin2m = np.sin(2.0 * math.pi / 180 * (m))

    Etime = (
            y * sin2l0
            - 2.0 * e * sinm
            + 4.0 * e * y * sinm * cos2l0
            - 0.5 * y * y * sin4l0
            - 1.25 * e * e * sin2m
    )

    return 180 / math.pi * Etime * 4.0


def eccentric_location_earth_orbit_vec(juliancentury: float) -> float:
    """Calculate the eccentricity of Earth's orbit"""
    return 0.016708634 - juliancentury * (0.000042037 + 0.0000001267 * juliancentury)


def var_y_vec(juliancentury: float) -> float:
    epsilon = obliquity_correction_vec(juliancentury)
    y = np.tan(math.pi / 180.0 * epsilon / 2.0)
    return y ** 2


def hour_angle_vec(
        latitude: np.array, declination: float, zenith: float, direction: SunDirection
) -> np.array:
    """Calculate the hour angle of the sun
    See https://en.wikipedia.org/wiki/Hour_angle#Solar_hour_angle
    Args:
        latitude: The latitude of the obersver
        declination: The declination of the sun
        zenith: The zenith angle of the sun
        direction: The direction of traversal of the sun
    Raises:
        ValueError
    """

    latitude_rad = latitude * math.pi / 180
    declination_rad = declination * math.pi / 180
    zenith_rad = zenith * math.pi / 180

    # n = cos(zenith_rad)
    # d = cos(latitude_rad) * cos(declination_rad)
    # t = tan(latitude_rad) * tan(declination_rad)
    # h = (n / d) - t

    h = (np.cos(zenith_rad) - np.sin(latitude_rad) * np.sin(declination_rad)) / (
            np.cos(latitude_rad) * np.cos(declination_rad)
    )

    HA = np.arccos(h)
    if direction == SunDirection.SETTING:
        HA = -HA
    return HA
