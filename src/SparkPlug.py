import numpy as np


class SparkPlug:

    """
    The SparkPlug class specifies properties of a spark plug where it is simply modeled
    the spark by a small region of temperature change
    """

    def __init__(self, pos, radius: float, temp: float, cycle: int, theta: float):

        # Initialize instance variables
        self.x = np.array(pos)
        self.radius = float(radius)
        self.temp = float(temp)
        self.cycle = int(cycle)
        self.theta = theta

    def aabb(self):
        """
        Determine the particles axis aligned bounding box
        :return: np.array of vertices
        """
        radius = self.radius
        return np.array([self.x[0] + radius, self.x[1] + radius],
                        [self.x[0] - radius, self.x[1] - radius])


pos = [1, 1]
radius = 1
temp = 2000
cycle = 4
theta = 3

derp = SparkPlug(pos, radius, temp, cycle, theta)
print(derp.cycle)
