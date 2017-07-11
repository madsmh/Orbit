import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import celestial3d
import read_horizon


# Physical properties of the celestial bodies in SI units
SUN_MASS = 1.98855e30
SUN_RADIUS = 695700000.0

EARTH_MASS = 5.97219e24
EARTH_RADIUS = 6371000.0

MERCURY_MASS = 3.302e23
MERCURY_RADIUS = 2439700

VENUS_MASS = 48.685e23
VENUS_RADIUS = 6051800.0

MARS_MASS = 6.4185e23
MARS_RADIUS = 3389500.0

JUPITER_MASS = 1898.13e24
JUPUTER_RADIUS = 71492000.0

SATURN_MASS = 5.68319e26
SATURN_RADIUS = 7149200

URANUS_MASS = 86.8103e24
URANUS_RADIUS = 24973000

NEPTUNE_MASS = 102.41e24
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
n = 64*400

# Time interval in seconds (1/64 Earth solar day)
dt = 1440*60/64


def check_mars():
    # Read the trajetory of Mars from file

    trajectory_mars = np.loadtxt("trajectories/mars.gz", float, delimiter=',')
    trajectory_earth = np.loadtxt("trajectories/earth.gz", float, delimiter=',')

    mars_nasa = np.array(mars_data)
    earth_nasa = np.array(earth_data)

    print(np.abs(trajectory_earth[:64:64*31, :] - earth_nasa[:, 0:3]))
    #print(np.abs(trajectory_mars[:len(mars_nasa), :]-mars_nasa[:, 0:3]))

def gen_coords():
    # Generate coordinates and save to file

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

    np.savetxt("trajectories/earth.gz", trajectory_earth, delimiter=',')
    np.savetxt("trajectories/mars.gz", trajectory_mars, delimiter=',')
    np.savetxt("trajectories/venus.gz", trajectory_venus, delimiter=',')
    np.savetxt("trajectories/jupiter.gz", trajectory_jupiter, delimiter=',')
    np.savetxt("trajectories/mercury.gz", trajectory_mercury, delimiter=',')
    np.savetxt("trajectories/sun.gz", trajectory_sun, delimiter=',')
    np.savetxt("trajectories/saturn.gz", trajectory_saturn, delimiter=',')
    np.savetxt("trajectories/uranus.gz", trajectory_uranus, delimiter=',')
    np.savetxt("trajectories/neptune.gz", trajectory_neptune, delimiter=',')
    np.savetxt("trajectories/pluto.gz", trajectory_pluto, delimiter=',')


def plot_planets():
    # Load coordinate data from files
    trajectory_earth = np.loadtxt("trajectories/earth.gz", float, delimiter=',')
    trajectory_mars = np.loadtxt("trajectories/mars.gz", float, delimiter=',')
    trajectory_venus = np.loadtxt("trajectories/venus.gz", float, delimiter=',')
    # trajectory_jupiter = np.loadtxt("trajectories/jupiter.gz", float, delimiter=',')
    trajectory_mercury = np.loadtxt("trajectories/mercury.gz", float, delimiter=',')
    trajectory_sun = np.loadtxt("trajectories/sun.gz", float, delimiter=',')
    # trajectory_saturn = np.loadtxt("trajectories/saturn.gz", float, delimiter=',')
    # trajectory_uranus = np.loadtxt("trajectories/uranus.gz", float, delimiter=',')
    # trajectory_neptune = np.loadtxt("trajectories/neptune.gz", float, delimiter=',')
    # trajectory_pluto = np.loadtxt("trajectories/pluto.gz", float, delimiter=',')

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # Plot inner planets
    ax.plot(trajectory_earth[:, 0], trajectory_earth[:, 1], trajectory_earth[:, 2], label=earth.name)
    ax.plot(trajectory_mars[:, 0], trajectory_mars[:, 1], trajectory_mars[:, 2], label=mars.name)
    ax.plot(trajectory_mercury[:, 0], trajectory_mercury[:, 1], trajectory_mercury[:, 2], label=mercury.name)
    ax.plot(trajectory_venus[:, 0], trajectory_venus[:, 1], trajectory_venus[:, 2], label=venus.name)

    # Sphere
    def plot_sphere(x0, y0, z0, body, col):
        """Plot a sphere at x0, y0, z0 representing body with color col"""

        u = np.linspace(0, 2 * np.pi, 20)
        v = np.linspace(0, np.pi, 20)
        x = x0 + body.radius * np.outer(np.cos(u), np.sin(v))
        y = y0 + body.radius * np.outer(np.sin(u), np.sin(v))
        z = z0 + body.radius * np.outer(np.ones(np.size(u)), np.cos(v))

        ax.plot_surface(x, y, z, color=col)

    # Plot the Sun
    plot_sphere(trajectory_sun[-1][0], trajectory_sun[-1][1], trajectory_sun[-1][2], sun, 'y')

    xylim = 1e11
    zlim = 1e11

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.auto_scale_xyz([-xylim, xylim], [-xylim, xylim], [-zlim, zlim])

    ax.legend()
    plt.show()
gen_coords()
check_mars()
plot_planets()
