import numpy as np
import matplotlib.pyplot as plt
import celestial

# Physical properties of the celestial bodies in SI units
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
MERCURY_MASS = 3.3011e23
MERCURY_RADIUS = 2439700

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

SATURN_X0 = -1352550000000.0
SATURN_Y0 = 0.0
SATURN_VX0 = 0.0
SATURN_VY0 = -10180.0
SATURN_MASS = 568.34e24
SATURN_RADIUS = 7149200

URANUS_X0 = -2741300000000.0
URANUS_Y0 = 0.0
URANUS_VX0 = 0.0
URANUS_VY0 = -7110
URANUS_MASS = 86.813e24
URANUS_RADIUS = 24973000

NEPTUNE_X0 = -4444450000000.0
NEPTUNE_Y0 = 0.0
NEPTUNE_VX0 = 0.0
NEPTUNE_VY0 = -5500
NEPTUNE_MASS = 102.413e24
NEPTUNE_RADIUS = 24341000

# Bodies
sun = celestial.Body('Sun', 0, 0, 0, 0, SUN_MASS, SUN_RADIUS)
earth = celestial.Body('Earth', EARTH_X0, EARTH_Y0, EARTH_VX0, EARTH_VY0, EARTH_MASS, EARTH_RADIUS)
jupiter = celestial.Body('Jupiter', JUPITER_X0, JUPITER_Y0, JUPITER_VX0, JUPITER_VY0, JUPITER_MASS, JUPUTER_RADIUS)
mars = celestial.Body('Mars', MARS_X0, MARS_Y0, MARS_VX0, MARS_VY0, MARS_MASS, MARS_RADIUS)
venus = celestial.Body('Venus', VENUS_X0, VENUS_Y0, VENUS_VX0, VENUS_VY0, VENUS_MASS, VENUS_RADIUS)
mercury = celestial.Body('Mercury', MERCURY_X0, MERCURY_Y0, MERCURY_VX0, MARS_VY0, MERCURY_MASS, MERCURY_RADIUS)
saturn = celestial.Body('Saturn', SATURN_X0, SATURN_Y0, SATURN_VX0, SATURN_VY0, SATURN_MASS, SATURN_RADIUS)
uranus = celestial.Body('Uranus', URANUS_X0, URANUS_Y0, URANUS_VX0, URANUS_VY0, URANUS_MASS, URANUS_RADIUS)
neptune = celestial.Body('Neptune', NEPTUNE_X0, NEPTUNE_Y0, NEPTUNE_VX0, NEPTUNE_VY0, NEPTUNE_MASS, NEPTUNE_RADIUS)

# Arrays of Bodies
sun_array = [mercury, earth, jupiter, mars, venus, saturn, uranus, neptune]
earth_array = [sun, jupiter, mars, venus, mercury, saturn, uranus, neptune]
jupiter_array = [sun, earth, mars, venus, mercury, saturn, uranus, neptune]
mars_array = [sun, earth, jupiter, venus, mercury, saturn, uranus, neptune]
venus_array = [sun, earth, jupiter, mars, mercury, saturn, uranus, neptune]
mercury_array = [sun, earth, jupiter, mars, venus, saturn, uranus, neptune]
saturn_array = [sun, earth, jupiter, mars, venus, mercury, uranus, neptune]
uranus_array = [sun, earth, jupiter, mars, venus, mercury, saturn, neptune]
neptune_array = [sun, earth, jupiter, mars, venus, mercury, saturn, uranus]

# Number of coordinate pairs
n = 8*600

# Time interval in seconds (1/8 Earth day)
dt = 86400.0/8

# NumPy arrays for coordinates
trajectory_earth = np.zeros(shape=(n, 2))
trajectory_mars = np.zeros(shape=(n, 2))
trajectory_venus = np.zeros(shape=(n, 2))
trajectory_jupiter = np.zeros(shape=(n, 2))
trajectory_mercury = np.zeros(shape=(n, 2))
trajectory_sun = np.zeros(shape=(n, 2))
trajectory_saturn = np.zeros(shape=(n, 2))
trajectory_uranus = np.zeros(shape=(n, 2))
trajectory_neptune = np.zeros(shape=(n, 2))


# Generate coordinates
def gen_coords():
    for i in range(n):
        trajectory_earth[i][0] = earth.x
        trajectory_earth[i][1] = earth.y

        trajectory_mars[i][0] = mars.x
        trajectory_mars[i][1] = mars.y

        trajectory_venus[i][0] = venus.x
        trajectory_venus[i][1] = venus.y

        trajectory_jupiter[i][0] = jupiter.x
        trajectory_jupiter[i][1] = jupiter.y

        trajectory_mercury[i][0] = mercury.x
        trajectory_mercury[i][1] = mercury.y

        trajectory_sun[i][0] = sun.x
        trajectory_sun[i][1] = sun.y

        trajectory_saturn[i][0] = saturn.x
        trajectory_saturn[i][1] = saturn.y

        trajectory_uranus[i][0] = uranus.x
        trajectory_uranus[i][1] = uranus.y

        trajectory_neptune[i][0] = neptune.x
        trajectory_neptune[i][1] = neptune.y

        earth.step(dt, earth_array)
        mars.step(dt, mars_array)
        venus.step(dt, venus_array)
        jupiter.step(dt, jupiter_array)
        mercury.step(dt, mercury_array)
        sun.step(dt, sun_array)
        saturn.step(dt, saturn_array)
        uranus.step(dt, uranus_array)
        neptune.step(dt, neptune_array)

gen_coords()

# Coordinates of the planets
x_earth, y_earth = trajectory_earth.T
x_mars, y_mars = trajectory_mars.T
x_venus, y_venus = trajectory_venus.T
x_jupiter, y_jupiter = trajectory_jupiter.T
x_mercury, y_mercury = trajectory_mercury.T
x_sun, y_sun = trajectory_sun.T
x_saturn, y_saturn = trajectory_saturn.T
x_uranus, y_uanus = trajectory_uranus.T
x_neptune, y_neptune = trajectory_neptune.T

plt.axes().set_aspect('equal', 'datalim')

plt.plot(x_earth, y_earth, 'g', linewidth=0.5)
plt.plot(x_mars, y_mars, 'r', linewidth=0.5)
plt.plot(x_venus, y_venus, 'b', linewidth=0.5)
plt.plot(x_jupiter, y_jupiter, 'k', linewidth=0.5)
plt.plot(x_mercury, y_mercury, 'k', linewidth=0.5)
plt.plot(x_sun, y_sun, 'k')
plt.plot(x_saturn, y_saturn, 'b', linewidth=0.5)
plt.plot(x_uranus, y_uanus, 'k', linewidth=0.5)
plt.plot(x_neptune, y_neptune, 'r', linewidth=0.5)

# Plot circles representing the planets
phi = np.linspace(0.0, 2*np.pi, 100)
na = np.newaxis

x_line_sun = x_sun[-1] + sun.radius * np.sin(phi[:, na])
y_line_sun = y_sun[-1] + sun.radius * np.cos(phi[:, na])

x_line_earth = x_earth[-1] + earth.radius * np.sin(phi[:, na])
y_line_earth = y_earth[-1] + earth.radius * np.cos(phi[:, na])

x_line_jupiter = x_jupiter[-1] + jupiter.radius * np.sin(phi[:, na])
y_line_jupiter = y_jupiter[-1] + jupiter.radius * np.cos(phi[:, na])


plt.plot(x_line_earth, y_line_earth)
plt.plot(x_line_sun, y_line_sun, 'y')
plt.plot(x_line_jupiter, y_line_jupiter, 'b')
plt.plot(x_line_sun, y_line_sun, 'y')


plt.show()

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
