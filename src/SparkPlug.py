import numpy as np


class SparkPlug:

    """
    The SparkPlug class specifies properties of a spark plug where it is simply modeled
    the spark by a small region of temperature change
    """

    def __init__(self, pos, radius, temp, cycle, theta):

        # Initialize instance variables
        self.pos = np.array(pos)
        self.radius = float(radius)
        self.temp = float(temp)
        self.cycle = cycle
        self.theta = theta


pos = [1, 1]
radius = 1
temp = 2000
cycle = 4
theta = 3

derp = SparkPlug(pos, radius, temp, cycle, theta)
print(derp.cycle)
