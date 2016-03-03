import numpy as np

class Wall:

    """
    Wall class is meant to represent the boundaries of the combustion chamber. Here walls
    are considered inelastic, with basic thermal properties. Wall boundaries are 2 dimensional
    represented by a vector of two vertices with 3 dimensional properties such as thickness
    passed as arguments.
    """

    def __init__(self, vert, rest_temp, conduct, cp, rho, d, area):

        # Initialize instance variables
        self.vert = vert
        self.rest_temp = float(rest_temp)
        self.temp = float(rest_temp)
        self.conduct = float(conduct)
        self.cp = float(cp)
        self.d = float(d)
        self.area = float(area)
        self.rho = float(rho)
        self.norm = self.calc_norm()

    def temp_change(self, q, t):

        # Heat loss rate assuming constant specific heat
        dt_heat = -self.conduct*self.area*(self.temp-self.rest_temp)/self.d + q
        self.temp += dt_heat*t/self.cp

    def calc_norm(self):

        # Calculate wall normal, later determine wall intersection
        norm = self.vert[1, :] - self.vert[0, :]
        norm[1] = -norm[1]
        norm /= np.linalg.norm(norm)
        return norm


class MovingWall(Wall):

    """
    The moving wall class inherits from the Wall class with added variables important for rigid
    body dynamics. Here Velocity is tracked and degrees of freedom are explicitly specified.
    """

    def __init__(self, vert, rest_temp, conduct, cp, rho, d, area, on_rails, tdof, rdof):

        # We're trying this inheritance thing
        Wall.__init__(self, vert, rest_temp, conduct, cp, rho, d, area)

        # Calculate derived values
        self.cen_mass = self.cen_mass()
        self.mass = self.calc_mass()
        self.connections = []
        self.vel = np.array([0, 0])
        self.on_rails = on_rails
        self.tdof = tdof
        self.rdof = rdof

    def calc_mass(self):
        # Determine mass of the wall
        length = np.linalg.norm(self.vert[1, :] - self.vert[0, :])
        mass = self.rho*self.area*self.d*length
        return mass

    def cen_mass(self):
        # Determine the center of mass
        cen_mass = (self.vert[0, :] + self.vert[1, :])/2
        return cen_mass


derp = np.array([[1, 1], [2, 2]])
herp = MovingWall(derp, 300, 1, 1, 1, 1, 1, 1, 1, 1)
merp = Wall(derp, 300, 1, 1, 1, 1, 1)
print(herp.mass)
herp.temp = 500
herp.temp_change(1, 1)
print(herp.temp)
print(herp.norm)
print(herp.cen_mass)
