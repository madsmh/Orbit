import numpy as np


class System:
    def __init__(self, bodies):

        self.bodies = bodies
        self.n = len(self.bodies)

    def get_positions(self):
        pos = np.zeros(self.n, 3)
        for b, i in zip(self.bodies, range(self.n)):
            pos[i][:] = b.get_position()
        return pos


class Body:
    """A celestial body class, with all initial values in SI units """

    def __init__(self, name, x0, y0, z0, vx0, vy0, vz0, mass, gm, radius):

        # Gravitational parameter
        self.GM = gm

        self.mass = mass

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

        # Radius of the body (m)
        self.radius = radius

    def compute_acceleration(self, x, y, z):
        """Computes the gravitational acceleration due to self at position (x, y) (m)"""
        # Deltas
        delta_x = self.x - x
        delta_y = self.y - y
        delta_z = self.z - z

        # Distance squared
        r2 = delta_x ** 2 + delta_y ** 2 + delta_z ** 2
        # Distance
        r = np.sqrt(r2)

        f = self.GM / r2

        # Return numpy array with the acceleration
        return np.array([f * delta_x / r, f * delta_y / r, f * delta_z / r])

    def get_position(self):
        return np.array([self.x, self.y, self.z])

    def ke(self):
        return 0.5 * self.mass * (self.vx ** 2 + self.vy ** 2 + self.vz ** 2)
