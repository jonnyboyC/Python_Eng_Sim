import numpy as np
import Quadtree


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
        self.aabb = Quadtree.Rectangle(
            self.x[0] - self.mole.radius,
            self.x[1] - self.mole.radius,
            2*self.mole.radius,
            2*self.mole.radius
        )
