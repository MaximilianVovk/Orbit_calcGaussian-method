"""
JD2GMST 
for a given Julian Date -> returns the angle between the First point of Aries
and Greenwich meridian (in radians)

XYZ2RADEC
for given x; y; z coordinates of a body in a geocentric equatorial coordinate
system with its origin at the Earth's centre, its z-axis along the north pole and its x-axis out
through the First Point of Aries -> returns the Right Ascension and Declination as seen
from the Earth's centre.

XYZ2CAR
for given x; y; z coordinates of a body in the geocentric equatorial coordinate
with the at the earth origin and for a given Julian date -> returns the Cartesian x; y; z coordinates in a geocentric
cartographic coordinate system with the same origin but with its z-axis along the north pole
and its x-axis out along the Greenwich meridian.

CAR2LATLON
accepts Cartesian x; y; z coordinates in the geocentric cartographic coordinate system
in km -> returns the geodetic latitude, longitude and altitude above the WGS84 reference
ellipsoid
"""

import numpy as np
import datetime

def JD2GMST(JD):
    """
    Computes the Greenwich Mean Sidereal Time (GMST) for a given Julian Date (JD).
    Input:
    JD: Julian Date
    Output:
    GMST: Greenwich Mean Sidereal Time (in radians)
    """
    # check if the input is a datetime object show an error message
    if isinstance(JD, datetime.datetime):
        raise TypeError('Input must be a Julian Date (float or int).')
    # check if the input is a string show an error message
    if isinstance(JD, str):
        raise TypeError('Input must be a Julian Date (float or int).')

    T = (JD - 2451545.0)/36525.0
    GMST = 280.46061837 + 360.98564736629*(T*36525.0) + 0.0003879332*T**2 - T**3/38710000
    GMST = GMST % 360
    GMST = GMST*np.pi/180.0
    return GMST

def XYZ2RADEC(x, y, z):
    """
    Computes the Right Ascension and Declination of a body in a geocentric equatorial
    coordinate system with its origin at the Earth's centre, its z-axis along the north
    pole and its x-axis out through the First Point of Aries.
    Input:
    x: x-coordinate of the body in the geocentric equatorial coordinate system
    y: y-coordinate of the body in the geocentric equatorial coordinate system
    z: z-coordinate of the body in the geocentric equatorial coordinate system
    GMST: Greenwich Mean Sidereal Time (in radians)
    Output:
    RA: Right Ascension (in radians)
    DEC: Declination (in radians)
    """
    RA = np.arctan2(y, x)*180.0/np.pi
    DEC = np.arcsin(z/np.sqrt(x**2 + y**2 + z**2))*180.0/np.pi
    return RA, DEC

def XYZ2CAR(x, y, z, GMST):
    """
    Computes the Cartesian x; y; z coordinates in a geocentric cartographic coordinate
    system with the same origin but with its z-axis along the north pole and its x-axis
    out along the Greenwich meridian.
    Input:
    x: x-coordinate of the body in the geocentric equatorial coordinate system
    y: y-coordinate of the body in the geocentric equatorial coordinate system
    z: z-coordinate of the body in the geocentric equatorial coordinate system
    GMST: Greenwich Mean Sidereal Time (in radians)
    Output:
    x: x-coordinate of the body in the geocentric cartographic coordinate system
    y: y-coordinate of the body in the geocentric cartographic coordinate system
    z: z-coordinate of the body in the geocentric cartographic coordinate system
    """
    x1 = x*np.cos(GMST) + y*np.sin(GMST)
    y1 = -x*np.sin(GMST) + y*np.cos(GMST)
    z1 = z
    return x1, y1, z1

def CAR2LATLON(x, y, z):
    """
    Computes the geodetic latitude, longitude and altitude above the WGS84 reference
    ellipsoid.
    Input:
    x: x-coordinate of the body in the geocentric cartographic coordinate system in km
    y: y-coordinate of the body in the geocentric cartographic coordinate system in km
    z: z-coordinate of the body in the geocentric cartographic coordinate system in km
    Output:
    lat: geodetic latitude (in radians)
    lon: longitude (in radians)
    alt: altitude above the WGS84 reference ellipsoid (in km)
    """
    # check if each commponents of the array magnitude is less than radius of the earth
    for i in range(len(x)):
        if np.sqrt(x[i]**2 + y[i]**2 + z[i]**2) < 6371:
            raise ValueError('Input must be in km and bigger than the radius of the Earth.')
        
    # transform in meters
    x=x*1000
    y=y*1000
    z=z*1000

    # WGS84 ellipsoid parameters
    a = 6378137.0
    b = 6356752.314245
    e2 = (a**2-b**2)/a**2
    e_2 = (a**2/b**2)-1
    p = np.sqrt(x**2 + y**2)
    u = np.arctan(z*a/(p*b))

    # compute latitude, longitude and altitude
    lat = np.arctan2((z + e_2*b*np.sin(u)**3),(p - e2*a*np.cos(u)**3))
    lon = np.arctan2(y, x)*180.0/np.pi
    alt = (p/np.cos(lat) - a/np.sqrt(1 - e2*np.sin(lat)**2))/1000
    lat = lat*180.0/np.pi
    return lat, lon, alt
