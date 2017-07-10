import numpy as np


class Body:
    """A celestial body class, with all initial values in SI units """

    def __init__(self, name, x0, y0, z0, vx0, vy0, vz0, mass, radius):

        # Constants of nature
        # Universal constant of gravitation
        self.G = 6.67408e-11

        # Name of the body (string)
        self.name = name

        # Initial position of the body (m)
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0

        # Position (m). Set to initial value.
        self.x = self.x0
        self.y = self.y0
        self.z = self.z0

        # Initial velocity of the body (m/s)
        self.vx0 = vx0
        self.vy0 = vy0
        self.vz0 = vz0

        # Velocity (m/s). Set to initial value.
        self.vx = self.vx0
        self.vy = self.vy0
        self.vz = self.vz0

        # Mass of the body (kg)
        self.M = mass

        # Radius of the body (m)
        self.radius = radius

    def compute_acceleration(self, x, y, z):
        """Computes the gravitational acceleration due to self at position (x, y) (m)"""

        # Acceleration in the x-direction (m/s^2)
        ax = self.G*self.M/((self.x-x)**2 + (self.y-y)**2 + (self.z-z)**2) * \
            (self.x-x)/np.sqrt((self.x-x)**2 + (self.y-y)**2 + (self.z-y)**2)

        # Acceleration in the y-direction (m/s^2)
        ay = self.G*self.M/((self.x-x)**2 + (self.y-y)**2 + (self.z-z)**2) * \
            (self.y-y)/np.sqrt((self.x-x)**2 + (self.y-y)**2 + (self.z-z)**2)

        # Acceleration in the z-direction (ms/s^2)
        az = self.G * self.M / ((self.x - x) ** 2 + (self.y - y) ** 2 + (self.z - z) ** 2) * \
            (self.z - z) / np.sqrt((self.x - x) ** 2 + (self.y - y) ** 2 + (self.z - z) ** 2)

        return ax, ay, az

    def step(self, dt, targets):
        """4th order Runge-Kutta integration"""

        # Acceleration due to targets (NumPy array)
        a = np.zeros(shape=(1, 3))
        for o in targets:
            a = a + np.array(o.compute_acceleration(self.x, self.y, self.z))

        k1x = self.vx
        k1y = self.vy
        k1z = self.vz

        k1vx = a[0][0]
        k1vy = a[0][1]
        k1vz = a[0][2]

        k2x = self.vx + dt / 2 * k1vx
        k2y = self.vy + dt / 2 * k1vy
        k2z = self.vz + dt / 2 * k1vz

        # Acceleration due to targets (NumPy array)
        a = np.zeros(shape=(1, 3))
        for o in targets:
            a += np.array(o.compute_acceleration(self.x + dt / 2 * k1x, self.y + dt / 2 * k1y,
                                                 self.z + dt / 2 * k1z))

        k2vx = a[0][0]
        k2vy = a[0][1]
        k2vz = a[0][2]

        k3x = self.vx + dt / 2 * k2vx
        k3y = self.vy + dt / 2 * k2vy
        k3z = self.vz + dt / 2 * k2vz

        # Acceleration due to targets (NumPy array)
        a = np.zeros(shape=(1, 3))
        for o in targets:
            a += np.array(o.compute_acceleration(self.x + dt / 2 * k2x, self.y + dt / 2 * k2y,
                                                 self.z + dt / 2 * k2z))

        k3vx = a[0][0]
        k3vy = a[0][1]
        k3vz = a[0][2]

        k4x = self.vx + dt * k3vx
        k4y = self.vy + dt * k3vy
        k4z = self.vz + dt * k3vz

        # Acceleration due to targets (NumPy array)
        a = np.zeros(shape=(1, 3))
        for o in targets:
            a += np.array(o.compute_acceleration(self.x + dt * k3x, self.y + dt * k3y,
                                                 self.z + dt * k3z))

        k4vx = a[0][0]
        k4vy = a[0][1]
        k4vz = a[0][2]

        # Update position
        self.x = self.x + dt / 6 * (k1x + 2 * k2x + 2 * k3x + k4x)
        self.y = self.y + dt / 6 * (k1y + 2 * k2y + 2 * k3y + k4y)
        self.z = self.z + dt / 6 * (k1z + 2 * k2z + 2 * k3z + k4z)

        # Update velocity
        self.vx = self.vx + dt / 6 * (k1vx + 2 * k2vx + 2 * k3vx + k4vx)
        self.vy = self.vy + dt / 6 * (k1vy + 2 * k2vy + 2 * k3vy + k4vy)
        self.vz = self.vz + dt / 6 * (k1vz + 2 * k2vz + 2 * k3vz + k4vz)
