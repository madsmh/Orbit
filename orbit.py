import numpy as np
# from scipy import integrate


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

        # Mass of the body (kg)
        self.M = mass

        # Radius of the body (m)
        self.radius = radius

    def compute_accelaraton(self, x, y):

        # Acceleration in the x-direction
        ax = self.G*self.M/((self.x-x)**2 + (self.y-y)**2)*(self.x-x)/np.sqrt((self.x-x)**2 + (self.y-y)**2)

        # Acceleration in the y-direction
        ay = self.G*self.M/((self.x-x)**2 + (self.y-y)**2)*(self.y-y)/np.sqrt((self.x-x)**2 + (self.y-y)**2)

        return ax, ay
