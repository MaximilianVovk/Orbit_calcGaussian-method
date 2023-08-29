"""
f2E(e,f)
return E
Converts true anomaly to eccentric anomaly.

E2f(e,E)
return f
Converts eccentric anomaly to true anomaly.

KEP2PER(µ, a, e, f)
return r_pqw, v_pqw
Converts Keplerian elements to perifocal coordinates.

PER2XYZ(r_pqw, v_pqw, i, Ω, ω)
return r, v
Converts perifocal coordinates to heliocentric Cartesian coordinates.

KEP2XYZ(µ, a, e, i, Ω, ω, f)
return r, v
Converts Keplerian elements to Cartesian coordinates.
"""

import numpy as np

def f2E(e,f):
    """
    Converts true anomaly to eccentric anomaly.

    Parameters:
        e (float): Eccentricity.
        f (float): True anomaly.

    Returns:
        E (float): Eccentric anomaly.
    """
    # if the orbit is parabolic or a hyperbolical raise an error
    if e>=1:
        raise ValueError('e eccentricity can be only e=0 or 0<e<1.')

    # Eccentric anomaly.
    E = 2 * np.arctan(np.sqrt((1 - e)/(1 + e))*np.tan(f/2))

    return E

def E2f(e,E):
    """
    Converts eccentric anomaly to true anomaly.

    Parameters:
        e (float): Eccentricity.
        E (float): Eccentric anomaly.

    Returns:
        f (float): True anomaly.
    """
    # if the orbit is parabolic or a hyperbolical raise an error
    if e>=1:
        raise ValueError('e eccentricity can be only e=0 or 0<e<1.')

    # True anomaly.
    f = 2 * np.arctan(np.sqrt((1 + e)/(1 - e))*np.tan(E/2))
    return f

def KEP2PER(µ, a, e, f):
    """
    Converts from Keplerian to perifocal coordinates.

    Parameters:
        µ (float): Gravitational parameter.
        a (float): Semi-major axis. (must be same unit as µ)
        e (float): Eccentricity.

    Returns:
        r (np.ndarray): Position vector.
        v (np.ndarray): Velocity vector.
    """
    # if µ is not defined correctly raise an error
    if µ<=0:
        raise ValueError('µ gravitational parameter must be µ>0.')
    # if the orbit is parabolic or a hyperbolical raise an error
    if e>=1:
        raise ValueError('e eccentricity can be only e=0 or 0<e<1.')

    # Semi-parameter.
    p = a * (1 - e**2)

    # Position in perifocal coordinates.
    r_pqw = np.array([p*np.cos(f)/(1 + e*np.cos(f)),
                      p*np.sin(f)/(1 + e*np.cos(f)),
                      0])

    # Velocity in perifocal coordinates.
    v_pqw = np.array([-np.sqrt(µ/p)*np.sin(f),
                      np.sqrt(µ/p)*(e + np.cos(f)),
                      0])

    return r_pqw, v_pqw


def PER2XYZ(r_pqw, v_pqw, i, Ω, ω):
    """
    Converts from perifocal to heliocentric Cartesian coordinates.

    Parameters:
        r_pqw (np.ndarray): Position vector in perifocal coordinates.
        v_pqw (np.ndarray): Velocity vector in perifocal coordinates.
        i (float): Inclination.
        Ω (float): Longitude of the ascending node.
        ω (float): Argument of periapsis.

    Returns:
        r (np.ndarray): Position vector in heliocentric coordinates.
        v (np.ndarray): Velocity vector in heliocentric coordinates.
    """

    # Transformation matrix from perifocal to inertial coordinates.
    C_pqw_to_i = np.array([[-np.sin(Ω)*np.cos(i)*np.sin(ω) + np.cos(Ω)*np.cos(ω),
                            -np.sin(Ω)*np.cos(i)*np.cos(ω) - np.cos(Ω)*np.sin(ω),
                            np.sin(Ω)*np.sin(i)],
                           [np.cos(Ω)*np.cos(i)*np.sin(ω) + np.sin(Ω)*np.cos(ω),
                            np.cos(Ω)*np.cos(i)*np.cos(ω) - np.sin(Ω)*np.sin(ω),
                            -np.cos(Ω)*np.sin(i)],
                           [np.sin(i)*np.sin(ω),
                            np.sin(i)*np.cos(ω),
                            np.cos(i)]])

    # Transformation to inertial coordinates.
    r = np.dot(C_pqw_to_i, r_pqw)
    v = np.dot(C_pqw_to_i, v_pqw)

    return r, v

def KEP2XYZ(µ, a, e, i, Ω, ω, f):
    """
    Converts from Keplerian to heliocentric Cartesian coordinates.

    Parameters:
        µ (float): Gravitational parameter.
        a (float): Semi-major axis. (must be same unit as µ)
        e (float): Eccentricity.
        i (float): Inclination.
        Ω (float): Longitude of the ascending node.
        ω (float): Argument of periapsis.
        f (float): True anomaly.

    Returns:
        r (np.ndarray): Position vector in heliocentric coordinates.
        v (np.ndarray): Velocity vector in heliocentric coordinates.
    """
    # Convert from Keplerian to perifocal coordinates.
    r_pqw, v_pqw = KEP2PER(µ, a, e, f)

    # Convert from perifocal to heliocentric Cartesian coordinates.
    r, v = PER2XYZ(r_pqw, v_pqw, i, Ω, ω)

    return r, v


