import numpy as np
import matplotlib.pyplot as plt

class Body:
    """A celestial body """

    def __init__(self, name, x0, y0, vx0, vy0, mass, radius):

        # Constants of nature
        self.G = 6.67408e-11

        # Name of the body (string)
        self.name = name

        # Initial position of the body (m)
        self.x0 = x0
        self.y0 = y0

        # Position (m)
        self.x = self.x0
        self.y = self.y0

        # Initial velocity of the body (m/s)
        self.vx0 = vx0
        self.vy0 = vy0

        # Velocity (m/s)
        self.vx = vx0
        self.vy = vy0

        # Mass of the body (kg)
        self.M = mass

        # Radius of the body (m)
        self.radius = radius

    def compute_acceleration(self, x, y):
        """Computes the gravitational acceleration due to self at position (x, y) (m)"""

        # Acceleration in the x-direction (m/s^2)
        ax = self.G*self.M/((self.x-x)**2 + (self.y-y)**2)*(self.x-x)/np.sqrt((self.x-x)**2 + (self.y-y)**2)

        # Acceleration in the y-direction (m/s^2)
        ay = self.G*self.M/((self.x-x)**2 + (self.y-y)**2)*(self.y-y)/np.sqrt((self.x-x)**2 + (self.y-y)**2)

        return ax, ay

    def step(self, dt, target):
        """4th order Runge-Kutta integration"""

        # Acceleration due to target (tuple)
        a = target.compute_acceleration(self.x, self.y)

        k1x = self.vx
        k1y = self.vy

        k1vx = a[0]
        k1vy = a[1]

        k2x = self.vx + dt / 2 * k1vx
        k2y = self.vy + dt / 2 * k1vy

        # Acceleration due to target (tuple)
        a = target.compute_acceleration(self.x + dt / 2 * k1x, self.y + dt / 2 * k1y)

        k2vx = a[0]
        k2vy = a[1]

        k3x = self.vx + dt / 2 * k2vx
        k3y = self.vy + dt / 2 * k2vy

        # Acceleration due to target (tuple)
        a = target.compute_acceleration(self.x + dt / 2 * k2x, self.y + dt / 2 * k2y)

        k3vx = a[0]
        k3vy = a[1]

        k4x = self.vx + dt * k3vx
        k4y = self.vy + dt * k3vy

        # Acceleration due to target (tuple)
        a = target.compute_acceleration(self.x + dt * k3x, self.y + dt * k3y)

        k4vx = a[0]
        k4vy = a[1]

        # Update position
        self.x = self.x + dt / 6 * (k1x + 2 * k2x + 2 * k3x + k4x)
        self.y = self.y + dt / 6 * (k1y + 2 * k2y + 2 * k3y + k4y)

        # Update velocity
        self.vx = self.vx + dt / 6 * (k1vx + 2 * k2vx + 2 * k3vx + k4vx)
        self.vy = self.vy + dt / 6 * (k1vy + 2 * k2vy + 2 * k3vy + k4vy)

SUN_MASS = 1.989e30
SUN_RADIUS = 695700000.0

EARTH_X0 = -147095000000.0
EARTH_Y0 = 0.0
EARTH_VX0 = 0.0
EARTH_VY0 = -30300.0
EARTH_MASS = 5.972e24
EARTH_RADIUS = 6371000.0
dt = 86400.0

# Bodies
sun = Body('Sun', 0, 0, 0, 0, SUN_MASS, SUN_RADIUS)
earth = Body('Earth', EARTH_X0, EARTH_Y0, EARTH_VX0, EARTH_VY0, EARTH_MASS, EARTH_RADIUS)

def plot_earth():
    """Generate x,y-coords of earth"""
    pos = np.zeros(shape=(365, 2))
    for i in range(365):
        pos[i][0] = earth.x
        pos[i][1] = earth.y
        earth.step(dt, sun)
    return pos

traject = plot_earth()
x, y = traject.T
r = earth.radius

phi = np.linspace(0.0, 2*np.pi, 100)

na = np.newaxis
x_line_earth = x[na, : ] + r * np.sin(phi[ : , na])
y_line_earth = y[na, : ] + r * np.cos(phi[ : , na])

x_line_sun = sun.radius*np.sin(phi[ : , na])
y_line_sun = sun.radius*np.cos(phi[ : , na])

plt.plot(x_line_earth, y_line_earth )
plt.plot(x_line_sun, y_line_sun)
plt.show()