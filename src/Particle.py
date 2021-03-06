import numpy as np
import random
import math
from Quadtree import Rectangle


class ChemicalReaction(object):
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
                    new_particles[i] = place(new_particles[0], new_particles[i - 1], new_particles[i])

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


class Molecule(object):
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


class Particle(object):
    """
    The Particle class specifies properties for a given particle, including which molecule
    it represents and its current position and velocity
    """
    g = np.array([0, -9.81])

    def __init__(self, molecule: Molecule, pos: list, vel: list):
        """
        Particle object represents an individual particles position and velocity
        :param molecule: Molecule object representing the chemical properties of the particle
        :param pos: np.array representing the position of the particle
        :param vel: np.array represting the velocity of the particle
        :return:
        """
        self.mole = molecule
        self.aabb = Rectangle(0, 0, 2*self.mole.radius, 2*self.mole.radius)
        self.x = pos
        self.u = np.array(vel)
        self.checked = False
        self.leaf = False

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, x):
        """
        Always update the bounding box when updating position
        :param x: update the internal _x and bounding box
        """
        self._x = np.array(x)
        self.aabb.x = x[0] - self.mole.radius,
        self.aabb.y = x[1] - self.mole.radius,


    def update_x(self, dt):
        """
        Update the particles position based on time step and velocity
        :param dt: time step
        :return: None
        """
        self.x += dt * self.u + 0.5 * Particle.g * dt ** 2

    def radiate_heat(self, dt):
        """
        Currently a stub
        """
        pass

    def __str__(self):
        """
        String representation of the particle
        :return:
        """
        return self.mole.symbol + " u mag = " + str(np.linalg.norm(self.u))



def place(p0: Particle, pi: Particle, pj: Particle):
    """
    Determines the center of particle of radius rj relative to tangent to particles p0 and pi. The new particle
    occupies the location clockwise relative to the vector defined centers of p0 and pi
    :rtype: Particle
    :param p0: Particle 0 at the center of the reaction cluster
    :param pi: Particle i surrounding and tangent to Particle 0
    :param rj: radius of Particle j
    :return: updated Particle j with new position
    """

    d = np.linalg.norm(p0.x - pi.x)
    rj = pj.mole.radius
    r0 = p0.mole.radius + rj
    ri = pi.mole.radius + rj
    x0 = p0.x
    xi = pi.x
    a = (r0 ** 2 - ri ** 2 + d ** 2) / (2 * d)
    h = np.sqrt(r0 ** 2 - a ** 2)
    p = x0 + np.dot(a / d, xi - x0)

    # New location
    pj.x[0] = p[0] + np.dot(h / d, xi[1] - x0[1])
    pj.x[1] = p[1] + np.dot(h / d, xi[0] - x0[0])
    return pj
