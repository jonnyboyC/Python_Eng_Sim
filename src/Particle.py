import numpy as np


class ChemReactions:

    """
    The chemical reactions class outlines the required reactants, products, activation energy, and chemical
    energy released or absorbed during the reaction
    """

    def __init__(self, reactants: list, products: list, stoch_r: list, stoch_p: list, e: float, de: float):
        self.react = reactants
        self.prod = products
        self.stoch_r = stoch_r
        self.stoch_p = stoch_p
        self.e = e
        self.de = de

    @property
    def reaction_set(self):
        # Return a set containing the symbols of the reactants

        reactants = set()
        for molecules in self.react:
            reactants.add(molecules.symbol)

        return reactants


class Molecule:

    """
    The molecule class, contains information on the elemental properties of a given molecule
    """

    def __init__(self, symbol: str, mass: float, radius: float):
        self.symbol = symbol
        self.mass = float(mass)
        self.radius = float(radius)


class Particle:
    """
    The particle class, specifies, properties for a given particle, including which molecule
    it represents and its current position and velocity
    """

    def __init__(self, molecule: Molecule, pos: list, vel: list):
        self.molecule = molecule
        self.pos = np.array(pos)
        self.vel = np.array(vel)

    def update_velocity(self, dt):
        # Update particles position
        self.pos += dt*self.vel

    def radiate_heat(self, dt):
        # Currently a stub, not sure if I'll implement this
        pass


O2 = Molecule("O_2", 5.3562e-26, 213e-12)
CH4 = Molecule("CH_4", 2.6781e-26, 234e-12)
CO = Molecule("CO", 4.6867e-26, 262e-12)
H2 = Molecule("H_2", 3.3476e-27, 157e-12)

fuel = ChemReactions([O2, CH4], [CO, H2], [1, 1], [1, 2], 100, 400)
print(fuel.reaction_set)
p1 = Particle(O2, [0, 0], [5, 1])
print(p1.pos)