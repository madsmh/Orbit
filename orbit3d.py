import numpy as np


class System:
    def __init__(self, properties):
        """Accepts a n x 10 array of initial properties for celestial bodies"""

        self.bodies = [Body(*p) for p in properties]
        self.n = len(self.bodies)

    def get_positions(self):
        """Returns a n x 3 array with coordinates"""
        pos = np.zeros(shape=(self.n, 3))
        for b, i in zip(self.bodies, range(self.n)):
            pos[i][:] = b.get_position()
        return pos

    def get_velocities(self):
        vel = np.zeros(shape=(self.n, 3))
        for b, i in zip(self.bodies, range(self.n)):
            vel[i][:] = b.get_velocity()
        return vel

    def get_accelerations(self):
        """Compute and return the resultant acceleration of
        all the bodies in an n x 3 array"""

        # Calculate the resultant force on each body
        resforce = np.sum(self.force_matrix(), axis=1)

        acc = np.zeros(shape=(self.n, 3))

        for a, i in zip(self.bodies, range(self.n)):
            acc[i][:] = a.compute_acceleration(resforce[i][:])

        return acc

    def set_positions(self, pos):
        """Accepts a n x 3 array with coordinates (x, y, z)"""
        for a, i in zip(self.bodies, range(self.n)):
            a.set_position(*pos[i][:])

    def set_velocities(self, vel):
        """Accepts a n x 3 array with velocities (vx, vy, vz)"""
        for a, i in zip(self.bodies, range(self.n)):
            a.set_position(*vel[i][:])


    def force_matrix(self):
        """Returns n x n x 3 array of all the forces in the system"""

        def force(body1, body2):
            """Vector force acting on body2 exerted by body1"""

            # Gravitational constant
            g = 6.67408e-11

            pos1 = body1.get_position()
            pos2 = body2.get_position()

            # Avoid a divide by zero
            if pos2 == pos1:
                return np.zeros(shape=3)

            r12 = pos2 - pos1
            dist = np.abs(r12)
            r12_hat = r12 / dist

            return -g * body1.mass * body2.mass / (dist ** 2) * r12_hat

        forces = np.zeros(shape=(self.n, self.n, 3))

        for one, i in zip(self.bodies, range(self.n)):
            for two, j in zip(self.bodies, range(self.n)):
                forces[i][j][:] = force(one, two)

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

        # Array to store integration results (x, y, x, vx, vy, vz)
        self.last_values = np.zeros(shape=(2, 6))

        # Integration counter
        self.i = 0

    def compute_acceleration(self, force):
        """Calculate the acceleration given the resultant force (vector)"""
        return force/self.mass

    def get_position(self):
        """Returns 1 x 3 array with the x, y, x positions"""
        return np.array([self.x, self.y, self.z])

    def get_velocity(self):
        """Returns a 1 x 3 array of velocities"""
        return np.array([self.vx, self.vy, self.vz])

    def get_last_values(self):
        """Returns the 2 x 3 array of the last two integration results"""
        return self.last_values

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
