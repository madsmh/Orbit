import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import celestial3d
import read_horizon
import os

# Physical properties of the celestial bodies in SI units
SUN_MASS = 1.98855e30
SUN_GM = 1.32712440018e20
SUN_RADIUS = 695700000.0

EARTH_MASS = 5.97219e24
EARTH_GM = 3.986004418e14
EARTH_RADIUS = 6371000.0

MERCURY_MASS = 3.302e23
MERCURY_GM = 22032.09 * 10 ** 9
MERCURY_RADIUS = 2439700

VENUS_MASS = 48.685e23
VENUS_GM = 324858.63 * 10 ** 9
VENUS_RADIUS = 6051800.0

MARS_MASS = 6.4185e23
MARS_GM = 4.282837e13
MARS_RADIUS = 3389500.0

JUPITER_MASS = 1898.13e24
JUPITER_GM = 126686511 * 10 ** 9
JUPUTER_RADIUS = 71492000.0

SATURN_MASS = 5.68319e26
SATURN_GM = 37931207.8 * 10 ** 9
SATURN_RADIUS = 7149200

URANUS_MASS = 86.8103e24
URANUS_GM = 5793966 * 10 ** 9
URANUS_RADIUS = 24973000

NEPTUNE_MASS = 102.41e24
NEPTUNE_GM = 6835107 * 10 ** 9
NEPTUNE_RADIUS = 24341000

PLUTO_MASS = 1.307e22
PLUTO_GM = 872.4 * 10 ** 9
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
sun = celestial3d.Body('Sun', 0, 0, 0, 0, 0, 0, SUN_GM, SUN_RADIUS)
earth = celestial3d.Body('Earth', *earth_data[0][:], EARTH_GM, EARTH_RADIUS)
jupiter = celestial3d.Body('Jupiter', *jupiter_data[0][:], JUPITER_GM, JUPUTER_RADIUS)
mars = celestial3d.Body('Mars', *mars_data[0][:], MARS_GM, MARS_RADIUS)
venus = celestial3d.Body('Venus', *venus_data[0][:], VENUS_GM, VENUS_RADIUS)
mercury = celestial3d.Body('Mercury', *mercury_data[0][:], MERCURY_GM, MERCURY_RADIUS)
saturn = celestial3d.Body('Saturn', *saturn_data[0][:], SATURN_GM, SATURN_RADIUS)
uranus = celestial3d.Body('Uranus', *uranus_data[0][:], URANUS_GM, URANUS_RADIUS)
neptune = celestial3d.Body('Neptune', *neptune_data[0][:], NEPTUNE_GM, NEPTUNE_RADIUS)
pluto = celestial3d.Body('Pluto', *pluto_data[0][:], PLUTO_GM, PLUTO_RADIUS)


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
detail = 1
n = detail * 700

# Time interval in seconds (1/detail Earth solar day)
dt = 86400 / detail


def check_mars():
    # Read the trajetory of Mars from file

    trajectory_mars = np.loadtxt("trajectories/mars.gz", float, delimiter=',')
    trajectory_earth = np.loadtxt("trajectories/earth.gz", float, delimiter=',')

    mars_nasa = np.array(mars_data)
    earth_nasa = np.array(earth_data)

    print('Absolule error from the NASA dataset (Earth): \n' +
          str(np.abs(trajectory_earth[0:31*detail:detail, :] - earth_nasa[:, 0:3])))
    print('Absolule error from the NASA dataset (Mars): \n' +
          str(np.abs(trajectory_mars[0:31*detail:detail, :] - mars_nasa[:, 0:3])))

    distance = (trajectory_earth[:, 0]-trajectory_mars[:, 0]) ** 2 + \
        (trajectory_earth[:, 1]-trajectory_mars[:, 1]) ** 2 + \
        (trajectory_earth[:, 2]-trajectory_mars[:, 2]) ** 2

    min_dist = np.sqrt(np.min(distance))
    print('Minimum distance between Earth and Mars is: ' + str(min_dist/1000) + ' km')

    orbit_dist = np.zeros(shape=2)
    dist = 0
    i = 0
    for x, y, z in trajectory_earth[detail * 300:detail * 400, :]:
        for u, v, w, in trajectory_mars[detail * 300:detail * 400, :]:
            if i == 0:
                orbit_dist[0] = (x-u) ** 2 + (y-v) ** 2 + (z-w) ** 2
                i += 1
            elif i == 1:
                orbit_dist[1] = (x-u) ** 2 + (y-v) ** 2 + (z-w) ** 2
                i += 1
            elif i > 1:
                dist = orbit_dist.min()
                orbit_dist[0] = dist
                orbit_dist[1] = (x - u) ** 2 + (y - v) ** 2 + (z - w) ** 2

    min_orbit_dist = np.sqrt(dist)

    print('Minium orbit distance is: ' + str(min_orbit_dist/1000) + ' km')


def gen_coords():
    if not os.path.exists('trajectories/'):
        os.makedirs('trajectories/')
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

    print('Trajectories generated and saved to file')


def plot_planets():
    # Load coordinate data from files
    trajectory_earth = np.loadtxt("trajectories/earth.gz", float, delimiter=',')
    trajectory_mars = np.loadtxt("trajectories/mars.gz", float, delimiter=',')
    trajectory_venus = np.loadtxt("trajectories/venus.gz", float, delimiter=',')
    trajectory_jupiter = np.loadtxt("trajectories/jupiter.gz", float, delimiter=',')
    trajectory_mercury = np.loadtxt("trajectories/mercury.gz", float, delimiter=',')
    trajectory_sun = np.loadtxt("trajectories/sun.gz", float, delimiter=',')
    trajectory_saturn = np.loadtxt("trajectories/saturn.gz", float, delimiter=',')
    trajectory_uranus = np.loadtxt("trajectories/uranus.gz", float, delimiter=',')
    trajectory_neptune = np.loadtxt("trajectories/neptune.gz", float, delimiter=',')
    trajectory_pluto = np.loadtxt("trajectories/pluto.gz", float, delimiter=',')

    mars_diag = np.array(read_horizon.readdiagnosticdata('mars'))
    earth_diag = np.array(read_horizon.readdiagnosticdata('earth'))

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # Plot planets
    #ax.plot(trajectory_mercury[:, 0], trajectory_mercury[:, 1], trajectory_mercury[:, 2], label=mercury.name)
    #ax.plot(trajectory_venus[:, 0], trajectory_venus[:, 1], trajectory_venus[:, 2], label=venus.name)
    ax.plot(trajectory_earth[:, 0], trajectory_earth[:, 1], trajectory_earth[:, 2], label='Earth simulated')
    ax.plot(trajectory_mars[:, 0], trajectory_mars[:, 1], trajectory_mars[:, 2], label='Mars simulated')

    ax.plot(earth_diag[:, 0], earth_diag[:, 1], earth_diag[:, 2], label='Earth diagnostic')
    ax.plot(mars_diag[:, 0], mars_diag[:, 1], mars_diag[:, 2], label='Mars diagnostic')

    # ax.plot(trajectory_jupiter[:, 0], trajectory_jupiter[:, 1], trajectory_jupiter[:, 2], label=jupiter.name)
    # ax.plot(trajectory_saturn[:, 0], trajectory_saturn[:, 1], trajectory_saturn[:, 2], label=saturn.name)
    # ax.plot(trajectory_uranus[:, 0], trajectory_uranus[:, 1], trajectory_uranus[:, 2], label=uranus.name)
    # ax.plot(trajectory_neptune[:, 0], trajectory_neptune[:, 1], trajectory_neptune[:, 2], label=neptune.name)
    # ax.plot(trajectory_pluto[:, 0], trajectory_pluto[:, 1], trajectory_pluto[:, 2], label=pluto.name)

    # Sphere
    def plot_sphere(x0, y0, z0, body, col):
        """Plot a sphere at x0, y0, z0 representing body with color col"""

        u = np.linspace(0, 2 * np.pi, 10)
        v = np.linspace(0, np.pi, 10)
        x = x0 + body.radius * np.outer(np.cos(u), np.sin(v))
        y = y0 + body.radius * np.outer(np.sin(u), np.sin(v))
        z = z0 + body.radius * np.outer(np.ones(np.size(u)), np.cos(v))

        ax.plot_surface(x, y, z, color=col)

    # Plot the Sun
    # plot_sphere(trajectory_sun[-1][0], trajectory_sun[-1][1], trajectory_sun[-1][2], sun, 'y')

    xylim = 1.2e11
    zlim = 1.2e11

    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.auto_scale_xyz([-xylim, xylim], [-xylim, xylim], [-zlim, zlim])

    ax.legend()
    plt.show()

gen_coords()
# check_mars()
plot_planets()
