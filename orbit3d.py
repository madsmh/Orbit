import numpy as np
# import matplotlib as mpl
# from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import celestial3d
import read_horizon


# Physical properties of the celestial bodies in SI units
SUN_MASS = 1.989e30
SUN_RADIUS = 695700000.0

EARTH_MASS = 5.972e24
EARTH_RADIUS = 6371000.0

MERCURY_MASS = 3.3011e23
MERCURY_RADIUS = 2439700

VENUS_MASS = 4.8675e24
VENUS_RADIUS = 6051800.0

MARS_MASS = 6.4185e23
MARS_RADIUS = 3389500.0

JUPITER_MASS = 1898.19e24
JUPUTER_RADIUS = 71492000.0

SATURN_MASS = 568.34e24
SATURN_RADIUS = 7149200

URANUS_MASS = 86.813e24
URANUS_RADIUS = 24973000

NEPTUNE_MASS = 102.413e24
NEPTUNE_RADIUS = 24341000

PLUTO_MASS = 1.307e22
PLUTO_RADIUS = 1195e3

earth_data = read_horizon.readdata('earth')
mercury_data = read_horizon.readdata('mercury')
venus_data = read_horizon.readdata('venus')
jupiter_data = read_horizon.readdata('jupiter')
saturn_data = read_horizon.readdata('saturn')
uranus_data = read_horizon.readdata('uranus')
neptune_data = read_horizon.readdata('neptune')
pluto_data = read_horizon.readdata('pluto')
mars_data = read_horizon.readdata('mars')

# Bodies
sun = celestial3d.Body('Sun', 0, 0, 0, 0, 0, 0, SUN_MASS, SUN_RADIUS)
earth = celestial3d.Body('Earth', *earth_data[0][:], EARTH_MASS, EARTH_RADIUS)
jupiter = celestial3d.Body('Jupiter', *jupiter_data[0][:], JUPITER_MASS, JUPUTER_RADIUS)
mars = celestial3d.Body('Mars', *mars_data[0][:], MARS_MASS, MARS_RADIUS)
venus = celestial3d.Body('Venus', *venus_data[0][:], VENUS_MASS, VENUS_RADIUS)
mercury = celestial3d.Body('Mercury', *mercury_data[0][:], MERCURY_MASS, MERCURY_RADIUS)
saturn = celestial3d.Body('Saturn', *saturn_data[0][:], SATURN_MASS, SATURN_RADIUS)
uranus = celestial3d.Body('Uranus', *uranus_data[0][:], URANUS_MASS, URANUS_RADIUS)
neptune = celestial3d.Body('Neptune', *neptune_data[0][:], NEPTUNE_MASS, NEPTUNE_RADIUS)
pluto = celestial3d.Body('Pluto', *pluto_data[0][:], PLUTO_MASS, PLUTO_RADIUS)


# Arrays of Bodies
sun_array = [mercury, earth, jupiter, mars, venus, saturn, uranus, neptune, pluto]
earth_array = [sun, jupiter, mars, venus, mercury, saturn, uranus, neptune, pluto]
jupiter_array = [sun, earth, mars, venus, mercury, saturn, uranus, neptune, pluto]
mars_array = [sun, earth, jupiter, venus, mercury, saturn, uranus, neptune, pluto]
venus_array = [sun, earth, jupiter, mars, mercury, saturn, uranus, neptune, pluto]
mercury_array = [sun, earth, jupiter, mars, venus, saturn, uranus, neptune, pluto]
saturn_array = [sun, earth, jupiter, mars, venus, mercury, uranus, neptune, pluto]
uranus_array = [sun, earth, jupiter, mars, venus, mercury, saturn, neptune, pluto]
neptune_array = [sun, earth, jupiter, mars, venus, mercury, saturn, uranus, pluto]
pluto_array = [sun, earth, jupiter, mars, venus, mercury, saturn, uranus, neptune]

# Number of coordinate pairs
n = 16*687

# Time interval in seconds (1/8 Earth day)
dt = 86400.0/16

# NumPy arrays for coordinates
trajectory_earth = np.zeros(shape=(n, 3))
trajectory_mars = np.zeros(shape=(n, 3))
trajectory_venus = np.zeros(shape=(n, 3))
trajectory_jupiter = np.zeros(shape=(n, 3))
trajectory_mercury = np.zeros(shape=(n, 3))
trajectory_sun = np.zeros(shape=(n, 3))
trajectory_saturn = np.zeros(shape=(n, 3))
trajectory_uranus = np.zeros(shape=(n, 3))
trajectory_neptune = np.zeros(shape=(n, 3))
trajectory_pluto = np.zeros(shape=(n, 3))


# Generate coordinates
def gen_coords():
    for i in range(n):
        trajectory_earth[i][0] = earth.x
        trajectory_earth[i][1] = earth.y
        trajectory_earth[i][2] = earth.z

        trajectory_mars[i][0] = mars.x
        trajectory_mars[i][1] = mars.y
        trajectory_mars[i][2] = mars.z

        trajectory_venus[i][0] = venus.x
        trajectory_venus[i][1] = venus.y
        trajectory_venus[i][2] = venus.z

        trajectory_jupiter[i][0] = jupiter.x
        trajectory_jupiter[i][1] = jupiter.y
        trajectory_jupiter[i][2] = jupiter.z

        trajectory_mercury[i][0] = mercury.x
        trajectory_mercury[i][1] = mercury.y
        trajectory_mercury[i][2] = mercury.z

        trajectory_sun[i][0] = sun.x
        trajectory_sun[i][1] = sun.y
        trajectory_sun[i][2] = sun.z

        trajectory_saturn[i][0] = saturn.x
        trajectory_saturn[i][1] = saturn.y
        trajectory_saturn[i][2] = saturn.z

        trajectory_uranus[i][0] = uranus.x
        trajectory_uranus[i][1] = uranus.y
        trajectory_uranus[i][2] = uranus.z

        trajectory_neptune[i][0] = neptune.x
        trajectory_neptune[i][1] = neptune.y
        trajectory_neptune[i][2] = neptune.z

        trajectory_pluto[i][0] = pluto.x
        trajectory_pluto[i][1] = pluto.y
        trajectory_pluto[i][2] = pluto.z

        earth.step(dt, earth_array)
        mars.step(dt, mars_array)
        venus.step(dt, venus_array)
        jupiter.step(dt, jupiter_array)
        mercury.step(dt, mercury_array)
        sun.step(dt, sun_array)
        saturn.step(dt, saturn_array)
        uranus.step(dt, uranus_array)
        neptune.step(dt, neptune_array)
        pluto.step(dt, pluto_array)

gen_coords()

# TODO Create x, y, z-coordinates for pluto
# Coordinates of the planets
x_earth, y_earth, z_earth = trajectory_earth.T
x_mars, y_mars, z_mars = trajectory_mars.T
x_venus, y_venus, z_venus = trajectory_venus.T
x_jupiter, y_jupiter, z_jupiter = trajectory_jupiter.T
x_mercury, y_mercury, z_mercury = trajectory_mercury.T
x_sun, y_sun, z_sun = trajectory_sun.T
x_saturn, y_saturn, z_saturn = trajectory_saturn.T
x_uranus, y_uanus, z_uranus = trajectory_uranus.T
x_neptune, y_neptune, z_neptune = trajectory_neptune.T

# TODO Write the 3D plotting for the planets and the sun.

fig = plt.figure()
ax = fig.gca(projection='3d')

ax.plot(x_earth, y_earth, z_earth, label=earth.name)
ax.plot(x_mars, y_mars, z_mars, label=mars.name)
# ax.plot(x_mercury, y_mercury, z_mercury, label=mercury.name)
ax.plot(x_venus, y_venus, z_venus, label=venus.name)
lim = 1.5e11

ax.auto_scale_xyz([-lim, lim], [-lim, lim], [-lim, lim])

ax.legend()
plt.show()

# plt.plot(x_mars, y_mars, 'r', linewidth=0.5)
# plt.plot(x_venus, y_venus, 'b', linewidth=0.5)
# plt.plot(x_jupiter, y_jupiter, 'k', linewidth=0.5)
# plt.plot(x_mercury, y_mercury, 'k', linewidth=0.5)
# plt.plot(x_sun, y_sun, 'k')
# plt.plot(x_saturn, y_saturn, 'b', linewidth=0.5)
# plt.plot(x_uranus, y_uanus, 'k', linewidth=0.5)
# plt.plot(x_neptune, y_neptune, 'r', linewidth=0.5)
#
# # Draw circles representing Mars
# x_line_mars = x_mars[na, :] + mars.radius * np.sin(phi[:, na])
# y_line_mars = y_mars[na, :] + mars.radius * np.cos(phi[:, na])
#
# # Draw circles representing Venus
# x_line_venus = x_venus[na, :] + venus.radius * np.sin(phi[:, na])
# y_line_venus = y_venus[na, :] + venus.radius * np.cos(phi[:, na])
#
# # Draw circles representing Venus
# x_line_jupiter = x_jupiter[na, :] + jupiter.radius * np.sin(phi[:, na])
# y_line_jupiter = y_jupiter[na, :] + jupiter.radius * np.cos(phi[:, na])
#
# # Draw circle representing the Sun
# x_line_sun = sun.radius * np.sin(phi[:, na])
# y_line_sun = sun.radius * np.cos(phi[:, na])
#
#
# plt.plot(x_line_earth, y_line_earth)
# plt.plot(x_line_sun, y_line_sun)
# plt.plot(x_line_mars, y_line_mars)
# plt.plot(x_line_venus, y_line_venus)
# plt.plot(x_line_jupiter, y_line_jupiter)
