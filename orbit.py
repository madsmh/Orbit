import numpy as np
import matplotlib.pyplot as plt
import celestial

# Physical properties of the celestial bodies
SUN_MASS = 1.989e30
SUN_RADIUS = 695700000.0

EARTH_X0 = -147095000000.0
EARTH_Y0 = 0.0
EARTH_VX0 = 0.0
EARTH_VY0 = -30300.0
EARTH_MASS = 5.972e24
EARTH_RADIUS = 6371000.0

MERCURY_X0 = -46.0e9
MERCURY_Y0 = 0.0
MERCURY_VX0 = 0.0
MERCURY_VY0 = -58980.0
MERCURY_MASS = 0.33011e24

VENUS_X0 = -107480000000.0
VENUS_Y0 = 0.0
VENUS_VX0 = 0.0
VENUS_VY0 = -35260.0
VENUS_MASS = 4.8675e24
VENUS_RADIUS = 6051800.0

MARS_X0 = -206620000000.0
MARS_Y0 = 0.0
MARS_VX0 = 0.0
MARS_VY0 = -26500.0
MARS_MASS = 6.4171e23
MARS_RADIUS = 3389500.0

JUPITER_X0 = -740520000000.0
JUPITER_Y0 = 0.0
JUPITER_VX0 = 0.0
JUPITER_VY0 = -13720.0
JUPITER_MASS = 1898.19e24
JUPUTER_RADIUS = 71492000.0

# Bodies
sun = celestial.Body('Sun', 0, 0, 0, 0, SUN_MASS, SUN_RADIUS)
earth = celestial.Body('Earth', EARTH_X0, EARTH_Y0, EARTH_VX0, EARTH_VY0, EARTH_MASS, EARTH_RADIUS)
jupiter = celestial.Body('Jupiter', JUPITER_X0, JUPITER_Y0, JUPITER_VX0, JUPITER_VY0, JUPITER_MASS, JUPUTER_RADIUS)
mars = celestial.Body('Mars', MARS_X0, MARS_Y0, MARS_VX0, MARS_VY0, MARS_MASS, MARS_RADIUS)
venus = celestial.Body('Venus', VENUS_X0, VENUS_Y0, VENUS_VX0, VENUS_VY0, VENUS_MASS, VENUS_RADIUS)

# Arrays of Bodies
earth_array = [sun, jupiter, mars, venus]
jupiter_array = [sun, earth, mars, venus]
mars_array = [sun, earth, jupiter, venus]
venus_array = [sun, earth, jupiter, mars]

# Number of coordinate pairs
n = 4333

# Time interval (1 Earth day)
dt = 86400.0

# NumPy arrays for coordinates
trajectory_earth = np.zeros(shape=(n, 2))
trajectory_mars = np.zeros(shape=(n, 2))
trajectory_venus = np.zeros(shape=(n, 2))
trajectory_jupiter = np.zeros(shape=(n, 2))

# Generate coordinates
for i in range(n):
    trajectory_earth[i][0] = earth.x
    trajectory_earth[i][1] = earth.y

    trajectory_mars[i][0] = mars.x
    trajectory_mars[i][1] = mars.y

    trajectory_venus[i][0] = venus.x
    trajectory_venus[i][1] = venus.y

    trajectory_jupiter[i][0] = jupiter.x
    trajectory_jupiter[i][1] = jupiter.y

    earth.step(dt, earth_array)
    mars.step(dt, mars_array)
    venus.step(dt, venus_array)
    jupiter.step(dt, jupiter_array)

# Coordinates of the planets
x_earth, y_earth = trajectory_earth.T
x_mars, y_mars = trajectory_mars.T
x_venus, y_venus = trajectory_venus.T
x_jupiter, y_jupiter = trajectory_jupiter.T

plt.axes().set_aspect('equal', 'datalim')

plt.plot(x_earth, y_earth, 'g', linewidth = 0.5)
plt.plot(x_mars, y_mars, 'r', linewidth = 0.5)
plt.plot(x_venus, y_venus, 'b', linewidth = 0.5)
plt.plot(x_jupiter, y_jupiter, 'r', linewidth = 0.5)
plt.show()

# phi = np.linspace(0.0, 2*np.pi, 100)
# na = np.newaxis
#
# # Draw circles representing Earth
# x_line_earth = x_earth[na, :] + earth.radius * np.sin(phi[:, na])
# y_line_earth = y_earth[na, :] + earth.radius * np.cos(phi[:, na])
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


