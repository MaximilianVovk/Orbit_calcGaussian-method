"""
position(JD, mu, a, e, i, OM, om, M,J0)
input JD and keplerian elements
returns the x,y,z coordinates of any planet in a heliocentric ecliptic frame.
by default the function returns the coordinates of the Earth in J2000.0

kepler(M, e)
input mean anomaly and eccentricity
returns the eccentric anomaly

diameter_from_magnitude(magnitude, distance_from_earth_vec, distance_from_sun_vec, albedo)
input the magnitude of the asteroid, the distance from the earth and the distance from the sun, and the albedo of the asteroid
returns the absolute magnitude and the diameter of the asteroid

rot_for_Earth_obliquity(vec, obliquity)
input a vector and the obliquity of the Planet
returns the vector rotated by the obliquity of the Planet, by default the obliquity of the Earth

gauss_method(obs)
input the observations of the object Julian Date, the Right Ascension in degrees and the Declination in degrees
returns the keplerian elements of the object and the covariance matrix
"""

import numpy as np
import math
from kep_xyz import E2f
from xyz_kep import XYZ2KEP
from kep_xyz import KEP2XYZ


def position(JD, mu=0.01720209895**2, a=1.00000011, e=0.01671022, i=np.deg2rad(0.00005),
              OM=np.deg2rad(-11.26064), om=np.deg2rad(102.94719+11.26064), 
              M=np.deg2rad(100.46435-(102.94719+11.26064)+11.26064),J0= 2451545.0):
    """
    input JD and keplerian elements 
    returns the x,y,z coordinates of any planet in a heliocentric ecliptic frame. 

    Parameters:
        JD (float): Julian date.
        mu (float): Gravitational parameter. (default mu in AU^3/day^2)
        a (float): Semi-major axis. (default Earth J2000.0)
        e (float): Eccentricity. (default Earth J2000.0)
        i (float): Inclination. (default Earth J2000.0)
        OM (float): Longitude of the ascending node. (default Earth J2000.0)
        om (float): Argument of perihelion. (default Earth J2000.0)
        M (float): Mean anomaly. (default Earth J2000.0)
        J0 (float): Julian date of the reference epoch. (default 2451545.0 = J2000.0)
    """
    # mean motion
    n = np.sqrt(mu/a**3)
    # calculate mean anomaly
    M = (M + n * (JD - J0))%( 2*np.pi )
    # calculate eccentric anomaly
    E = kepler(M, e)
    # calculate true anomaly
    f = E2f(e,E)
    # calculate x,y,z coordinates
    r, v = KEP2XYZ(mu, a, e, i, OM, om, f)
    return r

def kepler(M, e):
    """
    A subroutine to convert the eccentricty e and the mean anomaly M into the eccentric
    anomaly E the routine solve Kepler's equation iteratively. 

    Parameters:
        M (float): Mean anomaly.
        e (float): Eccentricity.

    Returns:
        E (float): Eccentric anomaly.
    """
    # calculate eccentric anomaly iterativelly
    E = M
    # a loop for 100 iterations
    for i in range(100):
        dE = (M - E+e*np.sin(E))/(1-e*np.cos(E))
        if abs(dE) < 1e-10:
            break
        E = E + dE
    return E

def diameter_from_magnitude(magnitude, distance_from_earth_vec, distance_from_sun_vec, albedo=0.1):
    """
    A subroutine that calculates the diameter of a body from its apparent magnitude, given its
    distance from Earth and the Sun and its albedo.

    Parameters:
        magnitude (float): Apparent magnitude of the body.
        distance_from_earth_vec (np.ndarray): Distance vector from the Earth to the body.
        distance_from_sun_vec (np.ndarray): Distance vector from the Sun to the body.
        albedo (float): Albedo of the body.

    Returns:
        diameter (float): Diameter of the body in km.
    """
    # Sun-object-Earth angle
    phase_ang = np.arccos(np.dot(distance_from_earth_vec, distance_from_sun_vec) /
                           (np.linalg.norm(distance_from_earth_vec) * np.linalg.norm(distance_from_sun_vec)))

    # asterid phase function
    phase_func = 2.5 * math.log10((1-0.15) * np.exp( -3.33*( np.tan(phase_ang/2))**0.63)+ 0.15*np.exp( -1.87*(np.tan(phase_ang/2))**1.22))

    # Convert magnitude to absolute magnitude for a photometric G parameter of 0.15.
    abs_magnitude = magnitude - 5 * math.log10(np.linalg.norm(distance_from_earth_vec) * 
                                               np.linalg.norm(distance_from_sun_vec)) + phase_func

    # Compute the diameter of the body
    diameter = 1329 * (10 ** (-abs_magnitude/5))/np.sqrt(albedo)

    return abs_magnitude, diameter

def magnitude(abs_magnitude, distance_from_earth_vec, distance_from_sun_vec):
    """
    A subroutine that calculates the apparent magnitude of a body from its diameter, given its absolute magnitude,
    distance from Earth and the Sun and its albedo.

    Parameters:
        abs_magnitude (float): Absolute magnitude of the body.
        distance_from_earth_vec (np.ndarray): Distance vector from the Earth to the body.
        distance_from_sun_vec (np.ndarray): Distance vector from the Sun to the body.

    Returns:
        magnitude (float): Apparent magnitude of the body.
    """
    # Sun-object-Earth angle
    phase_ang = np.arccos(np.dot(distance_from_earth_vec, distance_from_sun_vec) / 
                          (np.linalg.norm(distance_from_earth_vec) * np.linalg.norm(distance_from_sun_vec)))

    # asterid phase function
    phase_func = 2.5 * math.log10((1-0.15)*np.exp( -3.33*( np.tan(phase_ang/2))**0.63)+ 0.15* np.exp( -1.87*(np.tan(phase_ang/2))**1.22))

    # Convert magnitude to absolute magnitude for a photometric G parameter of 0.15.
    magnitude = abs_magnitude + 5 * math.log10(np.linalg.norm(distance_from_earth_vec) * 
                                               np.linalg.norm(distance_from_sun_vec)) - phase_func

    return magnitude

def rot_for_Earth_obliquity(vec, obliquity=-23.439291111111):
    """
    A subroutine that rotates a vector by the Earth's obliquity.

    Parameters:
        vec (np.ndarray): Vector to be rotated.
        obliquity (float): Earth's obliquity. (default Earth J2000.0)

    Returns:
        vec (np.ndarray): Rotated vector.
    """
    # Earth's obliquity
    obliquity = obliquity * np.pi/180

    # Rotation matrix
    R = np.array([[1, 0, 0], [0, np.cos(obliquity), np.sin(obliquity)], [0, -np.sin(obliquity), np.cos(obliquity)]])
    vec = np.dot(R, vec)

    return vec

def gauss_method(obs):
    """
    Converts from observations to orbital elements using Gauss' method using the Cunningham frame.

    Parameters:
        obs (np.ndarray): Julian Date, the Right Ascension in degrees and the Declination in degrees

    Returns:
        a (float): mean Semi-major axis. (in AU)
        e (float): mean Eccentricity.
        i (float): mean Inclination.
        Ω (float): mean Longitude of the ascending node.
        ω (float): mean Argument of periapsis.
        f (float): mean True anomaly.
        cov_kep (np.ndarray): Covariance matrix or Variance-covariance matrix of the orbital elements.

    """
    yr = 365.2568983263281 #Obs2OrbEle
    obliq= -23.439291111111 # Earth's obliquity
    mu = (4*np.pi**2)*(1+3.0404327497692654e-06) # SUN-Earth gravitational parameter in AU^3/y^2

    # orbiting body direction cosine vector
    rho_1= [np.cos(np.radians(obs[0][1]))*np.cos(np.radians(obs[0][2])), 
            np.sin(np.radians(obs[0][1]))*np.cos(np.radians(obs[0][2])), 
            np.sin(np.radians(obs[0][2]))]
    rho_2= [np.cos(np.radians(obs[1][1]))*np.cos(np.radians(obs[1][2])), 
            np.sin(np.radians(obs[1][1]))*np.cos(np.radians(obs[1][2])), 
            np.sin(np.radians(obs[1][2]))]
    rho_3= [np.cos(np.radians(obs[2][1]))*np.cos(np.radians(obs[2][2])), 
            np.sin(np.radians(obs[2][1]))*np.cos(np.radians(obs[2][2])), 
            np.sin(np.radians(obs[2][2]))]

    # Cunningham frame unit vectors
    Cun_1 = rho_1
    Cun_2 = np.cross(rho_1,np.cross(rho_3, rho_1)) / np.linalg.norm(np.cross(rho_1,np.cross(rho_3, rho_1)))
    Cun_3 = np.cross(Cun_1, Cun_2)

    Rot_Cun = np.array([[Cun_1[0],Cun_1[1],Cun_1[2]],
                        [Cun_2[0],Cun_2[1],Cun_2[2]],
                        [Cun_3[0],Cun_3[1],Cun_3[2]]
                        ])
    
    # Calculate the new vector in the Cunningham frame
    rho_1_Cun = np.dot(Rot_Cun,rho_1)
    rho_2_Cun = np.dot(Rot_Cun,rho_2)
    rho_3_Cun = np.dot(Rot_Cun,rho_3)

    # value of nu for the second observation
    if np.abs(rho_2_Cun[2]) < 1e-5:
        print('nu_2 value is close to 0 the calculation may be ill-conditioned.')

    # Define time intervals
    t_1 = (obs[2][0] - obs[1][0]) # 3-2
    t_2 = (obs[2][0] - obs[0][0]) # 3-1
    t_3 = (obs[1][0] - obs[0][0]) # 1-2

    # Define the sector ratios
    a1 = t_1/t_2 # it is b1
    a3 = t_3/t_2 # it is b3

    # Earth-Sun position vector in the heliocentric frame
    R_Earth_Helio_1=position(obs[0][0])
    R_Earth_1=rot_for_Earth_obliquity(R_Earth_Helio_1*(-1),obliq) 
    R_Earth_Helio_2=position(obs[1][0])
    R_Earth_2=rot_for_Earth_obliquity(R_Earth_Helio_2*(-1),obliq)
    R_Earth_Helio_3=position(obs[2][0])
    R_Earth_3=rot_for_Earth_obliquity(R_Earth_Helio_3*(-1),obliq)
    
    # Earth-Sun position vector in the cunningahm frame
    R_Earth_1=np.dot(Rot_Cun,R_Earth_1)
    R_Earth_2=np.dot(Rot_Cun,R_Earth_2)
    R_Earth_3=np.dot(Rot_Cun,R_Earth_3)

    # Calculate the distance of the object from the observer
    r_2 = (-a1*R_Earth_1[2]+R_Earth_2[2]-a3*R_Earth_3[2])/rho_2_Cun[2]
    r_3 = (r_2*rho_2_Cun[1]+a1*R_Earth_1[1]-R_Earth_2[1]+a3*R_Earth_3[1])/(a3*rho_3_Cun[1])
    r_1 = (r_2*rho_2_Cun[0]-a3*r_3*rho_3_Cun[0]+a1*R_Earth_1[0]-R_Earth_2[0]+a3*R_Earth_3[0])/a1

    # position in the in the cunningam frame
    R_1 = rot_for_Earth_obliquity(np.inner(Rot_Cun.T,np.array(
        [r_1*rho_1_Cun[0],r_1*rho_1_Cun[1],r_1*rho_1_Cun[2]])) - np.inner(Rot_Cun.T,R_Earth_1),-obliq)
    R_2 = rot_for_Earth_obliquity(np.inner(Rot_Cun.T,np.array(
        [r_2*rho_2_Cun[0],r_2*rho_2_Cun[1],r_2*rho_2_Cun[2]])) - np.inner(Rot_Cun.T,R_Earth_2),-obliq)
    R_3 = rot_for_Earth_obliquity(np.inner(Rot_Cun.T,np.array(
        [r_3*rho_3_Cun[0],r_3*rho_3_Cun[1],r_3*rho_3_Cun[2]])) - np.inner(Rot_Cun.T,R_Earth_3),-obliq)
    
    # Calculate the velocity of the object
    sigma_cost_1= mu/np.linalg.norm(R_1)**3
    sigma_cost_2= mu/np.linalg.norm(R_2)**3
    sigma_cost_3= mu/np.linalg.norm(R_3)**3

    # f parameters of Lagrange coefficients 
    f32 = 1 - 0.5*sigma_cost_2*(t_1/yr)**2
    f13 = 1 - 0.5*sigma_cost_3*(-t_2/yr)**2
    f21 = 1 - 0.5*sigma_cost_1*(t_3/yr)**2

    # g parameters of Lagrange coefficients
    g32 = t_1/yr - 1/6*sigma_cost_2*(t_1/yr)**3
    g13 = -t_2/yr - 1/6*sigma_cost_3*(-t_2/yr)**3
    g21 = t_3/yr - 1/6*sigma_cost_1*(t_3/yr)**3

    # velocity in the in the cunningam frame
    V_2 = (R_3-np.inner(f32,R_2))/g32
    V_3 = (R_1-np.inner(f13,R_3))/g13
    V_1 = (R_2-np.inner(f21,R_1))/g21

    # calculate the covarince matrix of the keplerian elements
    a_1, e_1, i_1, OM_1, om_1, f_1 = XYZ2KEP(mu,R_1,V_1)
    a_2, e_2, i_2, OM_2, om_2, f_2 = XYZ2KEP(mu,R_2,V_2)
    a_3, e_3, i_3, OM_3, om_3, f_3 = XYZ2KEP(mu,R_3,V_3)

    # Calculate the varance and covariance of x y and z and the velocity components
    a123=np.array([a_1,a_2,a_3])
    e123=np.array([e_1,e_2,e_3])
    i123=np.array([np.rad2deg(i_1),np.rad2deg(i_2),np.rad2deg(i_3)])
    OM123=np.array([np.rad2deg(OM_1),np.rad2deg(OM_2),np.rad2deg(OM_3)])
    om123=np.array([np.rad2deg(om_1),np.rad2deg(om_2),np.rad2deg(om_3)])
    f123=np.array([np.rad2deg(f_1),np.rad2deg(f_2),np.rad2deg(f_3)])

    data = np.array([a123, e123, i123,OM123,om123,f123])

    # Calculate the covariance matrix or variance-covariance matrix of the data
    cov_kep = np.cov(data, bias=True)

    return np.mean(a123), np.mean(e123), np.deg2rad(np.mean(i123)), np.deg2rad(np.mean(OM123)), np.deg2rad(np.mean(om123)), np.deg2rad(np.mean(f123)), cov_kep


