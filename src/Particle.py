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

    def activation(self, particles: list):
        # Determine if the relative motion of the particles produce enough kinetic energy to active the reaction
        sum_mass = sum(particles.mole.mass)
        sum_momentum = np.array([0, 0])
        for mass, vel in zip(particles.mole.mass, particles.vel):
            sum_momentum += mass*vel

        cen_vel = sum_momentum/sum_mass
        v0 = particles[0].vel - cen_vel
        v1 = particles[1].vel - cen_vel

        kinetic_energy =  1/2*sum_mass*np.linalg.norm(sum_momentum)**2
        kinetic_energy += 1/2*particles[0].molecule.mass*np.linalg.norm(v0)**2
        kinetic_energy += 1/2*particles[1].molecule.mass*np.linalg.norm(v1)**2

        if self.e > kinetic_energy:
            return particles




    def print(self):
        # Debugging function to show if reaction was read in correctly
        output = ""
        for molecule, stoch in zip(self.react, self.stoch_r):
            if stoch == 1:
                output += molecule.symbol + " + "
            else:
                output += str(stoch) + molecule.symbol + " + "

        # Remove trailing +
        output = output[:-3]
        output += " -> "

        for molecule, stoch in zip(self.prod, self.stoch_p):
            if stoch == 1:
                output += molecule.symbol + " + "
            else:
                output += str(stoch) + molecule.symbol + " + "

        # Remove trailing + and print
        output = output[:-3]
        print(output)


class Molecule:

    """
    The molecule class, contains information on the elemental properties of a given molecule
    """

    def __init__(self, symbol: str, mass: float, radius: float, inert: bool):
        self.symbol = symbol
        self.mass = float(mass)
        self.radius = float(radius)
        self.inert = bool(inert)


class Particle:
    """
    The particle class, specifies, properties for a given particle, including which molecule
    it represents and its current position and velocity
    """

    def __init__(self, molecule: Molecule, pos: list, vel: list):
        self.mole = molecule
        self.pos = np.array(pos)
        self.vel = np.array(vel)

    def update_velocity(self, dt):
        # Update particles position
        self.pos += dt*self.vel

    def radiate_heat(self, dt):
        # Currently a stub, not sure if I'll implement this
        pass

    def aabb(self):
        # Return Values for aspect adjusted bounding box.
        radius = self.mole.radius
        return np.array([self.pos[0] + radius, self.pos[1] + radius],
                        [self.pos[0] - radius, self.pos[1] - radius])


O2 = Molecule("O_2", 5.3562e-26, 213e-12)
CH4 = Molecule("CH_4", 2.6781e-26, 234e-12)
CO = Molecule("CO", 4.6867e-26, 262e-12)
H2 = Molecule("H_2", 3.3476e-27, 157e-12)

fuel = ChemReactions([O2, CH4], [CO, H2], [1, 1], [1, 2], 100, 400)
print(fuel.reaction_set)
p1 = Particle(O2, [0, 0], [5, 1])
print(p1.pos)
fuel.print()