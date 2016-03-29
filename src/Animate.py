from Combustion import step_time
from Quadtree import Quadtree
import numpy as np
import matplotlib


class Animate:

    def __init__(self, fig, walls):
        self.fig = fig
        self.fig.set_size_inches(16.0, 10.0)
        self.sim_ax = self.fig.add_axes([0.05, 0.05, 0.6, 0.9])
        self.temp_ax = self.fig.add_axes([0.7, 0.55, 0.25, 0.4])
        self.torq_ax = self.fig.add_axes([0.7, 0.05, 0.25, 0.4])
        self.line, = self.sim_ax.plot([], [], 'bo')
        self.w_line = []
        self.add_walls(walls)
        self.expand_bounds()

    def add_walls(self, walls: list):
        for wall in walls:
            line = self.sim_ax.plot(wall.vert[:,0], wall.vert[:,1], color='black')
            self.w_line.append(line)

    def init(self):
        self.line.set_data([],[])
        return self.line,

    def expand_bounds(self):
        limits = self.sim_ax.xaxis.get_data_interval()
        self.sim_ax.xaxis.set_view_interval(1.5*limits[0], 1.5*limits[1])

        limits = self.sim_ax.yaxis.get_data_interval()
        self.sim_ax.yaxis.set_view_interval(1.5*limits[0], 1.5*limits[1])

    def __call__(self, i, particles: list, objects: list, quad: Quadtree, dt: float):
        x = np.array([particle.x for particle in particles])
        self.line.set_data([x[:, 0], x[:, 1]])
        step_time(particles, objects, quad, dt)
        return self.line,