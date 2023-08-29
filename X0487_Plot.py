import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import seaborn as sn
import matplotlib.pyplot as plt

from kep_xyz import *
from xyz_kep import *
from radec_xyz import *
from Orbit_coordinates import *

mu = 4*np.pi**2 # SUN gravitational parameter in AU^3/y^2
def plot_Planet_orbit_3D_or_2D(Dim_3D_or_2D,n_planet=9):
    planet_names = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']
    planet = np.array([
        [0.38709893, 0.20563069, 7.00487, 48.33167, 77.45645, 252.25084],
        [0.72333199, 0.00677323, 3.39471, 76.68069, 131.53298, 181.97973],
        [1.00000011, 0.01671022, 0.00005, -11.26064, 102.94719, 100.46435],
        [1.52366231, 0.09341233, 1.85061, 49.57854, 336.04084, 355.45332],
        [5.20336301, 0.04839266, 1.30530, 100.55615, 14.75385, 34.40438],
        [9.53707032, 0.05415060, 2.48446, 113.71504, 92.43194, 49.94432],
        [19.19126393, 0.04716771, 0.76986, 74.22988, 170.96424, 313.23218],
        [30.06896348, 0.00858587, 1.76917, 131.72169, 44.97135, 304.88003],
        [39.48168677, 0.24880766, 17.14175, 110.30347, 224.06676, 238.92881],
    ])

    # plot the solar system north pole
    theta = np.pi/2
    x = 0
    y = 0
    z = 40
    # add a zero to the end of the array
    x = np.append(x, 0.0)
    y = np.append(y, 0.0)
    z = np.append(z, 0.0)
    ax.plot(x, y, z, 'k', markersize=10)

    # plot the first point of Aries as a line from the center
    x = 40
    y = 0
    z = 0
    # add a zero to the end of the array
    x = np.append(x, 0.0)
    y = np.append(y, 0.0)
    z = np.append(z, 0.0)
    ax.plot(x, y, z, 'k', markersize=10)

    ax.set_xlabel('y (AU)')
    ax.set_ylabel('x (AU)')
    ax.set_zlabel('z (AU)')

    ii=0
    # pik each argument in the each row of the n_planet
    for row in planet[0:n_planet]:
        a = row[0] # Semi-major axis in AU
        e = row[1] # Eccentricity 
        i = np.deg2rad(row[2]) # Inclination in radians
        OM = np.deg2rad(row[3]) # Longitude of the ascending node in radians
        om = np.deg2rad(row[4])-OM # Argument of periapsis in radians
        M = np.deg2rad(row[5])-OM-om # Mean anomaly in radians

        x_e=[]
        y_e=[]
        z_e=[]
        # plot the orbit with f from 1 to 360 deg
        for f in np.linspace(0, 360, 360):
            f = np.deg2rad(f) # True anomaly in radians
            r, v = KEP2XYZ(mu, a, e, i, OM, om, f)
            x_e = np.append(x_e, r[0])
            y_e = np.append(y_e, r[1])
            z_e = np.append(z_e, r[2])
            ax.set_zlabel('')

        # plot the orbit of the planet with the label of the planet_name
        if Dim_3D_or_2D==2:
            ax.plot(x_e, y_e,z_e, label=planet_names[ii])
            # view from the north pole
            ax.view_init(90, -90)
            # delete the z axis
            ax.set_zticks([])

        else:
            ax.plot(x_e, y_e, z_e, label=planet_names[ii])
        ii+=1
    # grid
    ax.grid(True)


# 1) ###################################################

# plot the orbits of the planets
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# 3D plot
plot_Planet_orbit_3D_or_2D(3)

Name_Ast='X0487'
obs = np.array([[2452465.5000000000000000,	    259.7657415651311794, 	-12.6153217571698129,  	18.047],
                [2452470.5000000000000000,   	259.2650934782980130, 	-12.6523752347416423, 	18.059],
                [2452480.5000000000000000,   	258.3754093732495676, 	-12.7520580832093593, 	18.085]])

# Name_Ast='X0506'
# obs = np.array([[2452465.5000000000000000,	1.3498908499056173,-17.0544601237825759, 34.232],
#                 [2452470.5000000000000000,   	1.3471518603560630,-17.0588726579291183, 34.231],
#                 [2452480.5000000000000000,   	1.3392532114351825,-17.0683461604096891, 34.229]])

# Name_Ast='X0613'
# obs = np.array([[2452465.5000000000000000,	    167.0760309907474266,1.1312156781805742, 23.377],
#                 [2452470.5000000000000000,   	167.1707119127673593,1.1081564268121049, 23.370],
#                 [2452480.5000000000000000,   	167.3941338167601032,1.0452617836226097, 23.353]])

# Name_Ast='X0718'
# obs = np.array([[2452465.5000000000000000, 15.6012629718496640, -9.4957892151581405, 25.010],
#                 [2452470.5000000000000000, 15.6120191227987597, -9.5276575325865682, 25.006],
#                 [2452480.5000000000000000, 15.6078224720019083, -9.6009907311125691, 24.996]])

# Name_Ast='X0769'
# obs = np.array([[2452465.5000000000000000,357.6833615090893659,0.3618335233054201, 11.951],
#                 [2452470.5000000000000000,357.8282418937238845,0.4280983288767462, 11.948],
#                 [2452480.5000000000000000,357.7820777378383355,0.4111137264269742, 11.928]])

a, e, i, OM, om, f, cov_kep = gauss_method(obs)

ax.set_title(Name_Ast+' || a='+str(np.round(a))+' e='+str(np.round(e,3))+' i='+str(np.round(np.rad2deg(i),1))+' OM='+str(np.round(np.rad2deg(OM),1))+' om='+str(np.round(np.rad2deg(om),1))+' f='+str(np.round(np.rad2deg(f),1)))
print('a=',a)
print('e=',e)
print('i=',np.rad2deg(i))
print('OM=',np.rad2deg(OM))
print('om=',np.rad2deg(om))
print('f=',np.rad2deg(f))
# # print('Error: ', cov_kep)

x_ast=[]
y_ast=[]
z_ast=[]
for f_i in np.linspace(0, 360, 360):
    f_i = np.deg2rad(f_i) # True anomaly in radians
    r, v = KEP2XYZ(mu, a, e, i, OM, om, f_i)
    x_ast = np.append(x_ast, r[0])
    y_ast = np.append(y_ast, r[1])
    z_ast = np.append(z_ast, r[2])

# plot the orbit of the planet with the label of the planet_name
ax.plot(x_ast, y_ast, z_ast, 'r', label=Name_Ast)

p_Ast,v_Ast = KEP2XYZ(mu, a, e, i, OM, om, f)
R_Earth_Helio_2=position(obs[1][0])

ax.plot(p_Ast[0], p_Ast[1], p_Ast[2], '.r', label=Name_Ast)
ax.plot(R_Earth_Helio_2[0], R_Earth_Helio_2[1], R_Earth_Helio_2[2], '.g', label='Earth JD')

# add a legend
ax.legend()
plt.show()

# plot the covariance matrix 
sn.heatmap(cov_kep, annot=True, fmt='g')
# label each tick
plt.xticks(np.arange(6)+0.5, ['a', 'e', 'i', 'OM', 'om', 'f'])
plt.yticks(np.arange(6)+0.5, ['a', 'e', 'i', 'OM', 'om', 'f'])
# show the plot

plt.show()

# 2) ###################################################

# plot the orbits of the planets
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# 3D plot
plot_Planet_orbit_3D_or_2D(2,6)

# print 2D plot of the orbit
ax.plot(x_ast, y_ast, 'r')
ax.plot(p_Ast[0], p_Ast[1], '.r', label=Name_Ast)

R_x_e=[]
R_y_e=[]
for ii in np.linspace(obs[1][0], obs[1][0]+365):
    R_Earth_current=position(ii)
    # put R_Earth_current in an array
    R_x_e=np.append(R_x_e, R_Earth_current[0])
    R_y_e=np.append(R_y_e, R_Earth_current[1])

ax.plot(R_x_e, R_y_e, 'g')
ax.plot(R_Earth_Helio_2[0], R_Earth_Helio_2[1], '.g', label='Earth')
# add a legend in the upper part of the plot and rotate it
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=3, fancybox=True, shadow=True)
ax.grid()
# dimensions of the x and y axis
ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)
ax.set_zlim(-5, 5)
plt.show()

distance_from_earth_vec = (p_Ast-R_Earth_Helio_2)
magnitude_2=obs[1][3]
albedo_min=0.02 #varies between 0.05 and 0.25,
abs_magnitude, diameter_min=diameter_from_magnitude(magnitude_2, distance_from_earth_vec, R_Earth_Helio_2, albedo_min)
albedo_MAX=0.06 #varies between 0.05 and 0.25,
abs_magnitude, diameter_MAX=diameter_from_magnitude(magnitude_2, distance_from_earth_vec, R_Earth_Helio_2, albedo_MAX)

print('Peri. distance = ',np.round(a*(1-e),2),'AU')
print('Diameter = ',np.round(diameter_min,2),'km - ',np.round(diameter_MAX,2),'km')
print('Abs. Magnitude = ',np.round(abs_magnitude,2))

# 3) ###################################################

# plot the orbits of the planets
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# 3D plot
plot_Planet_orbit_3D_or_2D(2,6)

E = f2E(e,f)

M = E - e*np.sin(E)

R_x_e2=[]
R_y_e2=[]
R_x_ast2=[]
R_y_ast2=[]
Mag=[]
RA=[]
DEC=[]
mu=0.01720209895**2
# print the position of the asteroid between 30 days and 120 days from obs[1][0]
for JDtime_Obs in np.linspace(obs[0][0]+30, obs[0][0]+120,91):

    R_Earth_Helio_all=position(JDtime_Obs)
    p_Ast_all=position(JDtime_Obs, mu, a, e, i, OM, om, M, obs[1][0])
    Mag=np.append(Mag, magnitude(abs_magnitude, p_Ast_all-R_Earth_Helio_all, R_Earth_Helio_all))
    R_Ast_from_Earth= rot_for_Earth_obliquity(p_Ast_all-R_Earth_Helio_all)
    
    # save position
    R_x_e2=np.append(R_x_e2, R_Earth_Helio_all[0])
    R_y_e2=np.append(R_y_e2, R_Earth_Helio_all[1])
    R_x_ast2=np.append(R_x_ast2, p_Ast_all[0])
    R_y_ast2=np.append(R_y_ast2, p_Ast_all[1])

    # convert the Cartesian coordinates to geodetic latitude, longitude and altitude
    RA_c, DEC_c = XYZ2RADEC(R_Ast_from_Earth[0], R_Ast_from_Earth[1], R_Ast_from_Earth[2])
    RA=np.append(RA, 360+RA_c)
    DEC=np.append(DEC, DEC_c)

# print 2D plot of the orbit
ax.plot(x_ast, y_ast, 'r')
ax.plot(R_x_ast2, R_y_ast2, '.r', label=Name_Ast)
ax.plot(R_x_e, R_y_e, 'g')
ax.plot(R_x_e2, R_y_e2, '.g', label='Earth')
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05), ncol=3, fancybox=True, shadow=True)
ax.grid()
ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)
ax.set_zlim(-5, 5)
plt.show()
plt.show()


# re create the time vector
JDtime_Obs=np.linspace(obs[0][0]+30, obs[0][0]+120,91)

# plot the magnitude and time of the asteroid
plt.plot(JDtime_Obs, Mag)
plt.xlabel('Time [JD]')
plt.ylabel('Magnitude')
plt.grid()
plt.show()

# getting the original colormap using cm.get_cmap() function
orig_map=plt.cm.get_cmap('viridis')
# reversing the original colormap using reversed() function
reversed_map = orig_map.reversed()
# plot the right ascension and declination of the asteroid with magnitude as color
plt.scatter(RA, DEC, c=Mag, cmap=reversed_map)
# invert the color bar
plt.colorbar()
plt.xlabel('Right Ascension [deg]')
plt.ylabel('Declination [deg]')
plt.grid()
plt.show()
