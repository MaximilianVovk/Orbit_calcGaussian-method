"""
the subrutines converts from cartesian to keplerian elements and viceversa

KEP2XYZ(µ, a, e, i, Ω, ω, f)
return r, v
Converts Keplerian elements to Cartesian coordinates.

XYZ2KEP(µ, r, v):
return a, e, i, Ω, ω, f
Converts Cartesian coordinates to Keplerian elements.
"""

import numpy as np



def KEP2XYZ(µ, a, e, i, Ω, ω, f):
    """
    Converts Keplerian elements to Cartesian coordinates.

    Parameters:
        µ (float): Gravitational parameter.
        a (float): Semi-major axis. (must be same unit as µ)
        e (float): Eccentricity.
        i (float): Inclination.
        Ω (float): Longitude of the ascending node.
        ω (float): Argument of periapsis.
        f (float): True anomaly.

    Returns:
        r (np.ndarray): Position vector.
        v (np.ndarray): Velocity vector.
    """
    # if µ is not defined correctly raise an error
    if µ<=0:
        raise ValueError('µ gravitational parameter must be µ>0.')
    # if is circular then the argument of periapsis is not defined.
    if e==0 and ω!=0:
        # raise a warning
        # raise ValueError('ω argument of periapsis is not defined for circular orbits, it must be set to zero ω=0.')
        print('WARNING ω argument of periapsis is not defined for circular orbits ω=0.')
        ω=0
    # if is equatorial then the longitude of the ascending node is not defined.
    if i==0 and Ω!=0:
        # raise a warning
        # raise ValueError('Ω longitude of the ascending node is not defined for equatorial orbits, it must be set to zero Ω=0.')
        print('WARNING Ω longitude of the ascending node is not defined for not inclined orbits Ω=0.')
        Ω=0
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

def XYZ2KEP(µ, r, v):
    """
    Converts Cartesian coordinates to Keplerian elements.

    Parameters:
        µ (float): Gravitational parameter.
        r (np.ndarray): Position vector. (must be same unit as µ)
        v (np.ndarray): Velocity vector. (must be same unit as µ)

    Returns:
        a (float): Semi-major axis.
        e (float): Eccentricity.
        i (float): Inclination.
        Ω (float): Longitude of the ascending node.
        ω (float): Argument of periapsis.
        f (float): True anomaly.
    """
    # if µ is not defined correctly raise an error
    if µ<=0:
        raise ValueError('µ gravitational parameter must be µ>0.')
    
    # Compute the specific angular momentum vector h
    h = np.cross(r, v)
    
    # Compute the magnitude of the specific angular momentum vector
    h_mag = np.linalg.norm(h)
    
    # Compute the unit vector of the specific angular momentum vector
    h_hat = h / h_mag
    
    # Compute the magnitude of the position vector r
    r_mag = np.linalg.norm(r)
    
    # Compute the velocity vector v
    v_mag = np.linalg.norm(v)
    
    # Compute the velocity vector v_hat
    v_hat = v / v_mag
    
    # Compute the eccentricity vector e_vec
    e_vec = (1/µ)*((v_mag**2 - µ/r_mag)*r - np.dot(r, v)*v)
    
    # Compute the eccentricity magnitude e
    e = np.linalg.norm(e_vec)
    
    # Compute the inclination i
    i = np.arccos(h_hat[2])

    # Compute the semi-major axis a
    a = 1 / (2/r_mag - v_mag**2/µ)

    # Compute the longitude of the ascending node Ω
    N = np.cross([0,0,1], h)
    N_mag = np.linalg.norm(N)

    if N_mag == 0:
        # Equatorial orbit
        Ω = 0.0
    else:
        # Non-equatorial orbit
        N_hat = N / N_mag
        Ω = np.arccos(N_hat[0])
        if N_hat[1] < 0:
            Ω = 2*np.pi - Ω
    
    # Compute the argument of periapsis ω    
    if e == 0:
        # Circular orbit
        ω = 0.0
    else:
        if i == 0:
            # Equatorial orbit
            ω = np.arctan2(e_vec[1], e_vec[0])
        else:
            # Non-circular orbit
            e_hat = e_vec / e
            ω = np.arccos(np.dot(N_hat, e_hat))
        if e_vec[2] < 0:
            ω = 2*np.pi - ω
    
    # Compute the true anomaly f
    if e == 0:
        if i == 0:
            # Circular equatorial orbit
            f = np.arccos(r[0] / r_mag)
        else:
            # Circular orbit
            f = np.arccos(np.dot(r, N) / (r_mag * N_mag))
    else:
        # Non-circular orbit
        cosE= (1/e)*(1 - r_mag/a)
        if cosE > 1:
            cosE = 1
        if cosE < -1:
            cosE = -1
        E = np.arccos(cosE)
        if np.dot(r, v) < 0:
            E = 2*np.pi - E
        f = 2*np.arctan(np.sqrt((1+e)/(1-e)) * np.tan(E/2))
        # f = np.arccos(np.dot(e_vec, r)/(e*r_mag))

    return a, e, i, Ω, ω, f
