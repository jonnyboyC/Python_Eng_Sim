import numpy as np
import random
import math


class ChemicalReaction:
    """
    The chemical reactions class outlines the required reactants, products, activation energy, and chemical
    energy released or absorbed during the reaction
    """

    def __init__(self, reactants: list, products: list, stoch_r: list, stoch_p: list, ea: float):
        """
        The ChemicalReaction Object represents the stoichiometric equation of a reaction,
        including the enthalpy of reaction and activation energy.
        :param reactants: List of Molecules constituting the reactants
        :param products: List of Molecules constituting the products
        :param stoch_r: List of Stoichiometric coefficients of the reactants
        :param stoch_p: List of Stoichiometric coefficients of the products
        :param ea: Float of activiation energy of the reaction
        :return:
        """
        self.react = reactants
        self.prod = products
        self.nu_r = stoch_r
        self.nu_p = stoch_p
        self.Ea = float(ea)
        self.Hr = self.delta_enthalpy()

    @property
    def reaction_set(self):
        # Return a set containing the symbols of the reactants

        reactants = set()
        for molecules in self.react:
            reactants.add(molecules.symbol)

        return reactants

    def delta_enthalpy(self):
        enthalpy_react = 0
        for react, nu in zip(self.react, self.nu_r):
            enthalpy_react += nu * react.Hf

        enthalpy_prod = 0
        for prod, nu in zip(self.prod, self.nu_p):
            enthalpy_prod += nu * prod.Hf

        return enthalpy_react - enthalpy_prod

    def activation(self, particles):
        # Determine if the relative motion of the particles produce enough kinetic energy to active the reaction

        # Determine the center of momentum frame
        sum_mass = sum(particle.mole.mass for particle in particles)
        sum_momentum = np.array([0.0, 0.0])
        sum_pos = np.array([0.0, 0.0])
        for particle in particles:
            sum_momentum += np.dot(particle.mole.mass, particle.u)
            sum_pos += np.dot(particle.mole.mass, particle.x)

        cen_mu = sum_momentum / sum_mass
        cen_m = sum_pos / sum_mass

        # Determine system kinetic energy in the center of momentum frame
        u = []
        for particle in particles:
            u.append(particle.u - cen_mu)

        kinetic_energy = 1 / 2 * sum_mass * np.linalg.norm(cen_mu) ** 2
        for particle, ui in zip(particles, u):
            kinetic_energy += 1 / 2 * particle.mole.mass * np.linalg.norm(ui) ** 2

        if self.Ea > kinetic_energy:
            # Need to decide if collision will occur here or up
            return particles
        else:
            total_energy = kinetic_energy + self.Hr
            new_particles = []
            for chem, nu in zip(self.prod, self.nu_p):
                for i in range(0, nu):
                    new_particles.append(Particle(chem, [0.0, 0.0], [0.0, 0, 0]))

            # Place new particles so center of mass is preserved
            random.shuffle(new_particles)
            if len(new_particles) > 1:
                theta = random.random() * 2 * math.pi
                sum_r = new_particles[0].radius + new_particles[1].radius
                new_particles[1].pos = np.dot(sum_r, [np.cos(theta), np.sin(theta)])
            if len(new_particles) > 2:
                for i in range(2, len(new_particles)):
                    new_particles[i].pos = place(new_particles[0], new_particles[i - 1], new_particles[i].radius)

            sum_mass = sum(particle.mole.mass for particle in particles)
            sum_pos = np.array([0.0, 0.0])
            for particle in new_particles:
                sum_pos += np.dot(particle.mole.mass, particle.x)

            del_cen_m = cen_m - sum_pos / sum_mass
            for particle in new_particles:
                particle.x += del_cen_m

            return new_particles

    def __str__(self):
        """
        Human readable representation of the reaction
        :return:
        """
        # Debugging function to show if reaction was read in correctly
        output = ""
        for molecule, stoch in zip(self.react, self.nu_r):
            if stoch == 1:
                output += molecule.symbol + " + "
            else:
                output += str(stoch) + molecule.symbol + " + "

        # Remove trailing +
        output = output[:-3]
        output += " -> "

        for molecule, stoch in zip(self.prod, self.nu_p):
            if stoch == 1:
                output += molecule.symbol + " + "
            else:
                output += str(stoch) + molecule.symbol + " + "

        # Remove trailing + and print
        output = output[:-3]
        return output

class Molecule:
    """
    The Molecule class, contains information on the elemental properties of a given molecule
    """
    avagondro = 6.0221409e+23

    def __init__(self, symbol: str, mass: float, radius: float, enthalpy: float, inert: bool):
        """
        The Molecule object represents the intrinsic properties of the particle
        :param symbol: String representing the chemical symbol
        :param mass:  Float representing the weight of the molecule in kg
        :param radius: Float representing a hard shell radius of the atom in m
        :param inert: Boolean representing if the gas is inert
        :param enthalpy: Float representing the heat of formation of the molecule
        :return:
        """
        self.symbol = str(symbol)
        self.mass = float(mass)
        self.radius = float(radius)
        self.Hf = self.per_particle(enthalpy)
        self.inert = bool(inert)

    @staticmethod
    def per_particle(enthalpy):
        return float(enthalpy) / (Molecule.avagondro * 1000)

    def __str__(self):
        return self.symbol


class Particle(Molecule):
    """
    The Particle class specifies properties for a given particle, including which molecule
    it represents and its current position and velocity
    """

    def __init__(self, molecule: Molecule, pos: list, vel: list):
        """
        Particle object represents an individual particles position and velocity
        :param molecule: Molecule object representing the chemical properties of the particle
        :param pos: np.array representing the position of the particle
        :param vel: np.array represting the velocity of the particle
        :return:
        """
        self.mole = molecule
        self.x = np.array(pos)
        self.u = np.array(vel)

    def update_velocity(self, dt):
        """
        Update the particles position based on time step and velocity
        :param dt: time step
        :return: None
        """
        self.x += dt * self.u

    def radiate_heat(self, dt):
        # Currently a stub, not sure if I'll implement this
        pass

    def aabb(self):
        """
        Determine the particles axis aligned bounding box
        :return: np.array of vertices
        """
        radius = self.radius
        return np.array([self.x[0] + radius, self.x[1] + radius],
                        [self.x[0] - radius, self.x[1] - radius])

    def __str__(self):
        """
        String representation of the particle
        :return:
        """
        return self.symbol + "u = " + np.linalg.norm(self.u)


def place(p0: Particle, pi: Particle, rj):
    """
    Determines the center of particle of radius rj relative to tangent to particles p0 and pi. The new particle
    occupies the location clockwise relative to the vector defined centers of p0 and pi
    :param p0: Particle 0 at the center of the reaction cluster
    :param pi: Particle i surrounding and tangent to Particle 0
    :param rj: radius of Particle j
    :return: np.array of the coordinates of the new center of Particle j
    """

    d = np.linalg.norm(p0.x - pi.x)
    r0 = p0.radius + rj
    ri = pi.radius + rj
    x0 = p0.x
    xi = pi.x
    a = (r0 ** 2 - ri ** 2 + d ** 2) / (2 * d)
    h = np.sqrt(r0 ** 2 - a ** 2)
    p = x0 + np.dot(a / d, xi - x0)

    # New location
    x = np.array([0.0, 0.0])
    x[0] = p[0] + np.dot(h / d, xi[1] - x0[1])
    x[1] = p[1] + np.dot(h / d, xi[0] - x0[0])
    return x

"""
O2 = Molecule("O_2", 5.3562e-26, 213e-12, 0, False)
CH4 = Molecule("CH_4", 2.6781e-26, 234e-12, -74.87, False)
CO = Molecule("CO", 4.6867e-26, 262e-12, -110.5, False)
H2 = Molecule("H_2", 3.3476e-27, 157e-12, 0, False)

fuel = ChemicalReaction([O2, CH4], [CO, H2], [1, 1], [1, 2], 1e-25)
print(fuel.reaction_set)
p1 = Particle(O2, [0, 0], [5, 1])
p2 = Particle(CO, [0, 0], [400, 1])
out = fuel.activation([p1, p2])
print(p1.x)
fuel.print()
"""