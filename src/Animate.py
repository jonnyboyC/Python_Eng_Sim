import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import numpy as np


class Animate:

    def __init__(self):
        self.fig = plt.figure(figsize=(16, 10))
        self.sim_ax = self.fig.add_axes([0.05, 0.05, 0.6, 0.9])
        self.temp_ax = self.fig.add_axes([0.7, 0.55, 0.25, 0.4])
        self.torq_ax = self.fig.add_axes([0.7, 0.05, 0.25, 0.4])
        self.w_line = []
        self.p_circle = []

    def show(self):
        self.fig.show()

    def add_walls(self, walls: list):
        for wall in walls:
            line = self.sim_ax.plot(wall.vert[:,0], wall.vert[:,1], color='black', linewidth=5)
            self.w_line.append(line)

    def add_particles(self, particles: list):
        for particle in particles:
            circle = Circle(particle.x, 1e8*particle.mole.radius)
            self.p_circle.append(circle)

        circles = PatchCollection(self.p_circle)
        self.sim_ax.add_collection(circles)
        pass

    def expand_bounds(self):
        limits = self.sim_ax.xaxis.get_data_interval()
        self.sim_ax.xaxis.set_view_interval(2*limits[0], 2*limits[1])

        limits = self.sim_ax.yaxis.get_data_interval()
        self.sim_ax.yaxis.set_view_interval(2*limits[0], 2*limits[1])




"""
fig, ax = generate_axis()
t = np.arange(0.0, 1.0, 0.01)
s = np.sin(2*np.pi*t)
line = ax.plot(t, s, color='blue', lw=2)
pass
"""