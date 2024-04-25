"""
This was pulled from https://github.com/kylebarron/suncalc-py/blob/master/suncalc/suncalc.py
And ported to MicroPython by Ethan Lacasse.

Port involved:
- Removing Numpy/Pandas functionality,
    - Pandas was not a dependancy, so removing Pandas
      mainly just meant deleting unused code.

    - Removing Numpy meant replacing all Numpy functions
      with their equivalent math functions,
      and replacing array math with list comprehension

- Adding epoch offset to convert between unix/embedded epoch standards

- Replacing soft Python constants with real MicroPython constants

- Making date parameters optional (uses time.now() if unset)

- Rewriting equations to use DecimalNumber for precision
  - This is much slower than the original math, however, MicroPython's single-precision floats
    gave practically useless results, only accurate to about half a day.
    This new approach should be obsurdly accurate (... and slow)
    and I think most embedded use cases for this module would prefer the accuracy.
  - mpy_decimal module has been forked and updated for use with this module

- Adding some additional docstrings and other formatting.


--------------------------


suncalc-py is ported from suncalc.js under the BSD-2-Clause license.

Copyright (c) 2014, Vladimir Agafonkin
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

# this is the difference between esp32 MicroPython epoch,
# and the unix epoch (2000-01-01 vs 1970-01-01)
# set to 0 if your platform uses unix epoch
_EPOCH_OFFSET = const(946684800)

# uses decimal module for high precision
# this is very slow, but very accurate.
# suncalc is not very useful with single-precision floats.
try:
    from .mpy_decimal import DecimalNumber as decimal
except:
    from suncalc.mpy_decimal import DecimalNumber as decimal



#from datetime import datetime
import time

#import numpy as np
import math

pd = None

# shortcuts for easier to read formulas
#PI = math.pi
PI = decimal.pi()
#sin = math.sin
#cos = math.cos
#tan = math.tan
#asin = math.asin
#atan = math.atan2
acos = math.acos
rad = PI / 180

def sin(decimal_input) -> decimal:
    if isinstance(decimal_input, (int, float)):
        return math.sin(decimal_input)
    return decimal_input.sin()

def cos(decimal_input) -> decimal:
    if isinstance(decimal_input, (int, float)):
        return math.cos(decimal_input)
    return decimal_input.cos()

def asin(decimal_input) -> decimal:
    if isinstance(decimal_input, (int, float)):
        return math.asin(decimal_input)
    return decimal_input.asin()

def tan(decimal_input) -> decimal:
    if isinstance(decimal_input, (int, float)):
        return math.tan(decimal_input)
    return decimal_input.tan()

def atan(input1, input2) -> decimal:
    
    if isinstance(input1, (int, float)):
        input1 = decimal(input1)
    if isinstance(input2, (int, float)):
        input2 = decimal(input2)
    
    return decimal.atan2(input1, input2)

def to_degrees(decimal_input) -> decimal:
    if isinstance(decimal_input, (int, float)):
        return math.degrees(decimal_input)
    return decimal_input.degrees()

def sign(in_num) -> int:
    if in_num == 0:
        return 0
    if in_num < 0:
        return -1
    if in_num >0:
        return 1

# def to_degrees(x):
#     if not isinstance(x, (int, float)):
#         x = float(x)
#     return math.degrees(x)
# sun times configuration (angle, morning name, evening name)
_DEFAULT_TIMES = (
    (-0.833, 'sunrise', 'sunset'),
    (-0.3, 'sunrise_end', 'sunset_start'),
    (-6, 'dawn', 'dusk'),
    (-12, 'nautical_dawn', 'nautical_dusk'),
    (-18, 'night_end', 'night'),
    (6, 'golden_hour_end', 'golden_hour'),
    )

# date/time constants and conversions
_DAY_MS = const(1000 * 60 * 60 * 24)
_DAY_SECONDS = const(60 * 60 * 24)
_J1970 = const(2440588)
_J2000 = const(2451545)


def to_milliseconds(epoch: int|float|tuple = None) -> decimal:
    """Get miliseconds from given epoch"""
    # original code returned time in milliseconds
    # using on datetime.timestamp()
    # For MicroPython, it's probably easier to just use the Time module
    
    # NOTE: this requires system time to be set to current UTC time!
    # Expected result is same as "datetime.now().timestamp()" in CPython (within 1000ms; int instead of float)
    
    # auto get time
    if epoch is None:
        epoch = time.time()
    
    # support for datetime tuples
    if type(epoch) == tuple:
        epoch = time.mktime(epoch)
    
    # adjust for unix vs embedded epoch
    epoch = decimal(epoch + _EPOCH_OFFSET)

    return epoch * 1000


def to_julian(date):
    return to_milliseconds(date) / _DAY_MS - 0.5 + _J1970


def from_julian(j):
    epoch = (j + 0.5 - _J1970) * _DAY_SECONDS
    
    # undo unix epoch offset
    epoch -= _EPOCH_OFFSET
    
    return time.localtime(int(epoch))


def to_days(date):
    return to_julian(date) - _J2000


# general calculations for position

# obliquity of the Earth
e = rad * decimal('23.4397')


def right_ascension(l, b):
    return atan(sin(l) * cos(e) - tan(b) * sin(e), cos(l))


def declination(l, b):
    return asin(sin(b) * cos(e) + cos(b) * sin(e) * sin(l))


def azimuth(H, phi, dec):
    return atan(sin(H), cos(H) * sin(phi) - tan(dec) * cos(phi))


def altitude(H, phi, dec):
    return asin(sin(phi) * sin(dec) + cos(phi) * cos(dec) * cos(H))


def sidereal_time(d, lw):
    return rad * (280.16 + (decimal('360.9856235') * d)) - lw


def astro_refraction(h):
    # the following formula works for positive altitudes only.
    # if h = -0.08901179 a div/0 would occur.
    if type(h) == list:
        h = [max(hi, 0) for hi in h]
        return [decimal('0.0002967') / tan(hi + decimal('0.00312536') / (hi + decimal('0.08901179'))) for hi in h]
    
    
    h = max(0, h)

    # formula 16.4 of "Astronomical Algorithms" 2nd edition by Jean Meeus
    # (Willmann-Bell, Richmond) 1998. 1.02 / tan(h + 10.26 / (h + 5.10)) h in
    # degrees, result in arc minutes -> converted to rad:
    return decimal('0.0002967') / tan(h + decimal('0.00312536') / (h + decimal('0.08901179')))


# general sun calculations


def solar_mean_anomaly(d):
    return rad * (decimal('357.5291') + decimal('0.98560028') * d)


def ecliptic_longitude(M):
    # equation of center
    C = rad * (decimal('1.9148') * sin(M) + 0.02 * sin(2 * M) + decimal('0.0003') * sin(3 * M))

    # perihelion of the Earth
    P = rad * decimal('102.9372')

    return M + C + P + PI


def sun_coords(d):
    M = solar_mean_anomaly(d)
    L = ecliptic_longitude(M)

    return {'dec': declination(L, 0), 'ra': right_ascension(L, 0)}


# calculations for sun times
J0 = decimal('0.0009')


def julian_cycle(d, lw):
    return round(d - J0 - lw / (2 * PI))


def approx_transit(Ht, lw, n):
    # no NP arrays, so we gotta accept regular lists
    if type(Ht) == list:
        return [J0 + (Hti + lw) / (2 * PI) + n for Hti in Ht]
    
    return J0 + (Ht + lw) / (2 * PI) + n


def solar_transit_j(ds, M, L):
    # no NP arrays, so we gotta accept regular lists
    if type(ds) == list:
        return [_J2000 + dsi + decimal('0.0053') * sin(M) - decimal('0.0069') * sin(2 * L) for dsi in ds]
    
    return _J2000 + ds + decimal('0.0053') * sin(M) - decimal('0.0069') * sin(2 * L)


def hour_angle(h, phi, d):
    # no NP arrays, so we gotta accept regular lists
    if type(h) == list:
        return [acos((sin(hi) - sin(phi) * sin(d)) / (cos(phi) * cos(d))) for hi in h]
    
    return acos((sin(h) - sin(phi) * sin(d)) / (cos(phi) * cos(d)))


def observer_angle(height):
    return decimal('-2.076') * math.sqrt(height) / 60


def get_set_j(h, lw, phi, dec, n, M, L):
    """Get set time for the given sun altitude
    """
    w = hour_angle(h, phi, dec)
    a = approx_transit(w, lw, n)
    return solar_transit_j(a, M, L)


def get_position(lng, lat, date=None, degrees=False):
    """Calculate sun position for a given date and latitude/longitude
    """
    
    lw = rad * -lng
    phi = rad * lat
    d = to_days(date)

    c = sun_coords(d)
    H = sidereal_time(d, lw) - c['ra']
    
    az = azimuth(H, phi, c['dec'])
    al = altitude(H, phi, c['dec'])
    
    
    if degrees:
        az = to_degrees(az)
        al = to_degrees(al)
    
    return {
        'azimuth': float(az),
        'altitude': float(al)}


def get_times(
        lng,
        lat,
        date=None,
        height=0,
        times: Iterable[Tuple[float, str, str]] = _DEFAULT_TIMES):
    """Calculate sun times

    Calculate sun times for a given date, latitude/longitude, and,
    optionally, the observer height (in meters) relative to the horizon
    """
    

    lw = rad * -lng
    phi = rad * lat

    dh = observer_angle(height)

    d = to_days(date)
    n = julian_cycle(d, lw)
    ds = approx_transit(0, lw, n)

    M = solar_mean_anomaly(ds)
    L = ecliptic_longitude(M)
    dec = declination(L, 0)

    Jnoon = solar_transit_j(ds, M, L)

    result = {
        'solar_noon': from_julian(Jnoon),
        'nadir': from_julian(Jnoon - 0.5)}

    angles = [time[0] for time in times]
    # h0 = (angles + dh) * rad
    h0 = [(angle + dh) * rad for angle in angles]
    
    # Need to add an axis for 2D broadcasting
    Jset = get_set_j(h0, lw, phi, dec, n, M, L)
    Jrise = [Jnoon - (Jseti - Jnoon) for Jseti in Jset]
    
    for idx, time in enumerate(times):
        result[time[1]] = from_julian(Jrise[idx])
        result[time[2]] = from_julian(Jset[idx])

    return result


# moon calculations, based on http://aa.quae.nl/en/reken/hemelpositie.html
# formulas


def moon_coords(d):
    """Geocentric ecliptic coordinates of the moon
    """

    # ecliptic longitude
    L = rad * (decimal('218.316') + decimal('13.176396') * d)
    # mean anomaly
    M = rad * (decimal('134.963') + decimal('13.064993') * d)
    # mean distance
    F = rad * (decimal('93.272') + decimal('13.229350') * d)

    # longitude
    l = L + rad * decimal('6.289') * sin(M)
    # latitude
    b = rad * decimal('5.128') * sin(F)
    # distance to the moon in km
    dt = 385001 - 20905 * cos(M)

    return {'ra': right_ascension(l, b), 'dec': declination(l, b), 'dist': dt}


def get_moon_position(lat, lng, date=None, degrees=False):

    lw = rad * -lng
    phi = rad * lat
    d = to_days(date)

    c = moon_coords(d)
    H = sidereal_time(d, lw) - c['ra']
    h = altitude(H, phi, c['dec'])

    # formula 14.1 of "Astronomical Algorithms" 2nd edition by Jean Meeus
    # (Willmann-Bell, Richmond) 1998.
    pa = atan(sin(H), tan(phi) * cos(c['dec']) - sin(c['dec']) * cos(H))

    # altitude correction for refraction
    h = h + astro_refraction(h)
    
    az = azimuth(H, phi, c['dec'])
    
    if degrees:
        az = to_degrees(az)
        h = to_degrees(h)

    return {
        'azimuth': az,
        'altitude': h,
        'distance': c['dist'],
        'parallacticAngle': pa}


# calculations for illumination parameters of the moon, based on
# http://idlastro.gsfc.nasa.gov/ftp/pro/astro/mphase.pro formulas and Chapter 48
# of "Astronomical Algorithms" 2nd edition by Jean Meeus (Willmann-Bell,
# Richmond) 1998.


def get_moon_illumination(date=None):

    d = to_days(date)
    s = sun_coords(d)
    m = moon_coords(d)

    # distance from Earth to Sun in km
    sdist = 149598000

    phi = acos(
        sin(s['dec']) * sin(m['dec']) +
        cos(s['dec']) * cos(m['dec']) * cos(s['ra'] - m['ra']))
    inc = atan(sdist * sin(phi), m['dist'] - sdist * cos(phi))
    angle = atan(
        cos(s['dec']) * sin(s['ra'] - m['ra']),
        sin(s['dec']) * cos(m['dec']) -
        cos(s['dec']) * sin(m['dec']) * cos(s['ra'] - m['ra']))

    return {
        'fraction': float((1 + cos(inc)) / 2),
        'phase': float(0.5 + 0.5 * inc * sign(angle) / PI),
        'angle': float(angle)}

