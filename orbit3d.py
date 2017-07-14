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

    def force_matrix(self):

        def force(body1, body2):
            """Force acting on body1 from body2"""

            g = 6.67408e-11
            pos1 = body1.get_position()
            pos2 = body2.get_position()

            # Avoid a divide by zero
            if pos1 == pos1:
                return np.array([0, 0, 0])

            r12 = pos2 - pos1
            dist = np.abs(r12)
            r12_hat = r12 / dist

            return -g * body1.mass * body2.mass / (dist ** 2) * r12_hat

        forces = np.zeros(shape=(self.n, self.n, 3))

        for one, i in zip(self.bodies, range(self.n)):
            for two, j in zip(self.bodies, range(self.n)):
                forces[j][i][:] = force(one, two)

        return forces


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

    def compute_acceleration(self, force):
        """Calculate the acceleration given the resultant force (vector)"""
        return force/self.mass

    def get_position(self):
        return np.array([self.x, self.y, self.z])

    def get_velocity(self):
        return np.array([self.vx, self.vy, self.vz])

    def set_position(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def set_velocity(self, vx, vy, vz):
        self.vx = vx
        self.vy = vy
        self.vz = vz

    def ke(self):
        """Calculate the kinetic energy"""
        return 0.5 * self.mass * (self.vx ** 2 + self.vy ** 2 + self.vz ** 2)
