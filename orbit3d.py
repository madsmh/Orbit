#     Solar system simulator. Copyright (c) 2017 Mads M. Hansen
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.


import numpy as np
import read_horizon
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import read_phys_properties as npp


class System:
    def __init__(self, names, posvel, gms, radii):
        """Accepts arrays of initial properties for celestial bodies"""

        self.n = len(names)

        # Initialize the bodies in our solar system
        self.bodies = [Body(names[i], *posvel[i], gms[i], radii[i]) for i in
                       range(self.n)]

    def get_positions(self):
        """Returns a n x 3 array with position coordinates"""
        pos = np.zeros(shape=(self.n, 3))
        for b, i in zip(self.bodies, range(self.n)):
            pos[i][:] = b.get_position()
        return pos

    def get_velocities(self):
        """Returns a n x 3 array with velocities"""
        vel = np.zeros(shape=(self.n, 3))
        for b, i in zip(self.bodies, range(self.n)):
            vel[i][:] = b.get_velocity()
        return vel

    def set_positions(self, pos):
        """Accepts a n x 3 array with coordinates (x, y, z)"""
        for a, i in zip(self.bodies, range(self.n)):
            a.set_position(*pos[i][:])

    def set_velocities(self, vel):
        """Accepts a n x 3 array with velocities (vx, vy, vz) and
        updates the positions of all bodies in this solar system"""
        for a, i in zip(self.bodies, range(self.n)):
            a.set_position(*vel[i][:])

    def get_accelerations(self):
        """Returns n x n array of the resultant accelerations in the system"""

        def acceleration(body1, body2):
            """Vector acceleration (in ms^-2) acting on body2 exerted by body1"""

            pos1 = body1.get_position()
            pos2 = body2.get_position()

            # Avoid a divide by zero
            if np.array_equal(pos2, pos1):
                return np.zeros(3)

            # Distance vector
            r12 = pos2 - pos1
            dist = np.linalg.norm(r12)

            return -body1.GM / (dist ** 3) * r12

        # Array to hold the force vectors
        accelerations = np.zeros(shape=(self.n, self.n, 3))

        # Calculate the force vectors and insert them into the array
        for a, i in zip(self.bodies, range(self.n)):
            for b, j in zip(self.bodies, range(self.n)):
                accelerations[i][j][:] = acceleration(a, b)

        accelerations = accelerations.sum(axis=0)

        return accelerations

class Body:
    """A celestial body class, with all initial values in SI units """


    def __init__(self, name, x0, y0, z0, vx0, vy0, vz0, gm, radius):

        # Gravitational parameter
        self.GM = gm

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

    def get_position(self):
        """Returns 1 x 3 array with the x, y, x positions"""
        return np.array([self.x, self.y, self.z])

    def get_velocity(self):
        """Returns a 1 x 3 array of velocities"""
        return np.array([self.vx, self.vy, self.vz])

    def set_position(self, x, y, z):
        """Set the position"""
        self.x = x
        self.y = y
        self.z = z

    def set_velocity(self, vx, vy, vz):
        """Set the velocity"""
        self.vx = vx
        self.vy = vy
        self.vz = vz


class Trajectory:
    """Saves 3-dimensional trajectories for a numbder of objects"""
    def __init__(self, n_trajectories, n_coords):
        self.trajectories = [np.zeros(shape=(n_coords, 3)) for _ in range(n_trajectories)]
        self.n_trajectories = n_trajectories
        self.row_counter = 0
        self.n_coords = n_coords

    def set_trajectory_position(self, pos):
        """Inputs a new position for every object from an n x 3 table"""
        for i in range(self.n_trajectories):
            self.trajectories[i][self.row_counter] = pos[i]
        self.row_counter += 1

    def get_trajectory(self, i):
        """Returns array i"""
        return self.trajectories[i]

    def get_position_at_index(self, i):
        """Gets the positions of all objects at index i
           as a n x 3 array"""

        data = np.zeros(shape=(self.n_trajectories, 3))

        for t, j in zip(self.trajectories, range(self.n_trajectories)):
            data[j] = t[i]

        return data
# List of body names
body_names, body_radii, body_gms = npp.read_phys_properties()
# (Mean-)Radii of the bodies (m)

n_bodies = len(body_names)

# Construct list of initial positions and velocities for each body (m and m/s)
init_pos_vel = np.zeros(shape=(n_bodies, 6))

for _, __ in zip(body_names, range(n_bodies)):
    init_pos_vel[__][:] = read_horizon.readdata(_.lower())[0]

# Solar system instance
detail = 1
dt = 86400/detail
n_rows = 1131*detail

sol = System(body_names, init_pos_vel, body_gms, body_radii)
tra = Trajectory(len(body_names), n_rows)

# Verlet


def verlet(system, trajectory, rows, delta_t):

    delta_t2 = delta_t ** 2

    # TODO Implement Velocity Verlet
    for k in range(rows):
        if k == 0:
            # Get initial positions
            q0 = system.get_positions()

            # Save to trajectory
            trajectory.set_trajectory_position(q0)

        elif k == 1:
            # Get previous position
            q0 = trajectory.get_position_at_index(0)

            # Get initial velocity
            p0 = system.get_velocities()

            # Calculate accerleration
            a = system.get_accelerations()

            # Calculate q1
            q1 = q0 + p0 * delta_t + 0.5 * a * delta_t2

            # Save to trajectory
            trajectory.set_trajectory_position(q1)

            # Update positions of the planets
            system.set_positions(q1)

        # Calculate q_n+1
        else:
            # Calculate accerleration
            a = system.get_accelerations()

            # Get the prevous results
            qn1 = tra.get_position_at_index(k-2)
            qn = tra.get_position_at_index(k-1)

            # Calculate new new positions
            qplus = 2*qn - qn1 + a * delta_t2

            # Save to trajectory
            trajectory.set_trajectory_position(qplus)

            # Update positions of the planets
            system.set_positions(qplus)


verlet(sol, tra, n_rows, dt)

# Plot the orbits

fig = plt.figure()
ax = fig.gca(projection='3d')

# venus = read_horizon.readdiagnosticdata('venus')
# luna_diagnostic = read_horizon.readdiagnosticdata('luna')


def diangnostic():
    earth_diagnostic = read_horizon.readdiagnosticdata('earth')[:, 0:3]
    earth_sim = tra.get_trajectory(3)
    sun_sim = tra.get_trajectory(0)

    # Coordinates of sim earth with respect to the sun
    earth_new_coords = (earth_sim - sun_sim)[::detail, :]

    error = np.linalg.norm(earth_new_coords/1000-earth_diagnostic/1000, axis=1)
    max_error = np.max(error)
    std = np.std(error)
    mean_error = np.mean(error)

    print('Mean error is: ' + str(mean_error) + ' km')
    print('Std. dev. is : ' + str(std))
    print('Max. error is: ' + str(max_error) + ' km')

    # print(list(zip(body_names, body_gms)))


# diangnostic()

for j in range(len(body_names)):
    ax.plot(tra.get_trajectory(j)[:, 0], tra.get_trajectory(j)[:, 1],
            tra.get_trajectory(j)[:, 2], label=body_names[j])

# ax.plot(venus[:, 0], venus[:, 1], venus[:, 2], label='Venus diagnostic')
# ax.plot(earth[:, 0], earth[:, 1], earth[:, 2], label='Earth diagnostic')
# ax.plot(luna_diagnostic[:, 0], luna_diagnostic[:, 1], luna_diagnostic[:, 2], label='Luna diagnostic')

dim = 1e12
ax.auto_scale_xyz([-dim, dim], [-dim, dim], [-dim, dim])
plt.legend()
plt.show()
