import numpy as np


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

    def step(self, dt, *targets):
        """4th order Runge-Kutta integration"""

        # Acceleration due to targets (NumPy array)
        a = 0
        for i in range(len(targets)):
            a += np.array(targets[i].compute_acceleration(self.x, self.y))

        k1x = self.vx
        k1y = self.vy

        k1vx = a[0]
        k1vy = a[1]

        k2x = self.vx + dt / 2 * k1vx
        k2y = self.vy + dt / 2 * k1vy

        # Acceleration due to targets (NumPy array)
        a = 0
        for i in range(len(targets)):
            a += np.array(targets[i].compute_acceleration(self.x + dt / 2 * k1x, self.y + dt / 2 * k1y))

        k2vx = a[0]
        k2vy = a[1]

        k3x = self.vx + dt / 2 * k2vx
        k3y = self.vy + dt / 2 * k2vy

        # Acceleration due to targets (NumPy array)
        a = 0
        for i in range(len(targets)):
            a += np.array(targets[i].compute_acceleration(self.x + dt / 2 * k2x, self.y + dt / 2 * k2y))

        k3vx = a[0]
        k3vy = a[1]

        k4x = self.vx + dt * k3vx
        k4y = self.vy + dt * k3vy

        # Acceleration due to targets (NumPy array)
        a = 0
        for i in range(len(targets)):
            a += np.array(targets[i].compute_acceleration(self.x + dt * k3x, self.y + dt * k3y))

        k4vx = a[0]
        k4vy = a[1]

        # Update position
        self.x = self.x + dt / 6 * (k1x + 2 * k2x + 2 * k3x + k4x)
        self.y = self.y + dt / 6 * (k1y + 2 * k2y + 2 * k3y + k4y)

        # Update velocity
        self.vx = self.vx + dt / 6 * (k1vx + 2 * k2vx + 2 * k3vx + k4vx)
        self.vy = self.vy + dt / 6 * (k1vy + 2 * k2vy + 2 * k3vy + k4vy)
