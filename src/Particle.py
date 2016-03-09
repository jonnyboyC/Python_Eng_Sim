import numpy as np
import random
import math


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

    def activation(self, particles):
        # Determine if the relative motion of the particles produce enough kinetic energy to active the reaction

        # Determine the center of momentum frame
        sum_mass = sum(particle.mole.mass for particle in particles)
        sum_momentum = np.array([0.0, 0.0])
        sum_pos = np.array([0.0, 0.0])
        for particle in particles:
            sum_momentum += np.dot(particle.mole.mass, particle.vel)
            sum_pos += np.dot(particle.mole.mass, particle.pos)

        cen_v = sum_momentum/sum_mass
        cen_m = sum_pos/sum_mass

        # Determine system kinetic energy in the center of momentum frame
        v = []
        for particle in particles:
            v.append(particle.vel - cen_v)

        kinetic_energy = 1/2*sum_mass*np.linalg.norm(cen_v)**2
        for particle, vi in zip(particles, v):
            kinetic_energy += 1/2*particle.mole.mass*np.linalg.norm(vi)**2

        if self.e > kinetic_energy:
            return particles
        else:
            total_energy = kinetic_energy + self.de
            new_particles = []
            for chem, stoch in zip(self.prod, self.stoch_p):
                for i in range(0, stoch):
                    new_particles.append(Particle(chem, [0.0,0.0], cen_v))

            # Place new particles so center of mass is preserved
            random.shuffle(new_particles)
            if len(new_particles) > 1:
                theta = random.random()*2*math.pi
                sum_r = new_particles[0].mole.radius + new_particles[1].mole.radius
                new_particles[1].pos = np.dot(sum_r, [np.cos(theta), np.sin(theta)])
            if len(new_particles) > 2:
                for i in range(2, len(new_particles)):
                    new_particles[i].pos = place(new_particles[0], new_particles[i-1], new_particles[i].mole.radius)

            sum_mass = sum(particle.mole.mass for particle in particles)
            sum_pos = np.array([0.0, 0.0])
            for particle in new_particles:
                sum_pos += np.dot(particle.mole.mass, particle.pos)

            del_cen_m = cen_m - sum_pos/sum_mass
            for particle in new_particles:
                particle.pos += del_cen_m

            return new_particles


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
    The Molecule class, contains information on the elemental properties of a given molecule
    """

    def __init__(self, symbol: str, mass: float, radius: float, inert: bool):
        self.symbol = symbol
        self.mass = float(mass)
        self.radius = float(radius)
        self.inert = bool(inert)


class Particle:
    """
    The Particle class specifies properties for a given particle, including which molecule
    it represents and its current position and velocity
    """

    def __init__(self, molecule: Molecule, pos: list, vel: list):
        self.mole = molecule
        self.pos = np.array(pos)
        self.vel = np.array(vel)

    def update_velocity(self, dt):
        """
        Update the particles position based on time step and velocity
        :param dt: time step
        :return: None
        """
        self.pos += dt*self.vel

    def radiate_heat(self, dt):
        # Currently a stub, not sure if I'll implement this
        pass

    def aabb(self):
        """
        Determine the particles axis aligned bounding box
        :return: np.array of vertices
        """
        radius = self.mole.radius
        return np.array([self.pos[0] + radius, self.pos[1] + radius],
                        [self.pos[0] - radius, self.pos[1] - radius])


def place(p0: Particle, pi: Particle, rj):
    """
    Determines the center of particle of radius rj relative to tangent to particles p0 and pi. The new particle
    occupies the location clockwise relative to the vector defined centers of p0 and pi
    :param p0: Particle 0 at the center of the reaction cluster
    :param pi: Particle i surrounding and tangent to Particle 0
    :param rj: radius of Particle j
    :return: np.array of the coordinates of the new center of Particle j
    """

    d = np.linalg.norm(p0.pos - pi.pos)
    r0 = p0.mole.radius+rj
    ri = p1.mole.radius+rj
    x0 = p0.pos
    xi = pi.pos
    a = (r0**2 - ri**2 + d**2)/(2*d)
    h = np.sqrt(r0**2 - a**2)
    p = x0 + np.dot(a/d, xi - x0)

    # New location
    pos = np.array([0.0,0.0])
    pos[0] = p[0] + np.dot(h/d, xi[1]-x0[1])
    pos[1] = p[1] + np.dot(h/d, xi[0]-x0[0])
    return pos


O2 = Molecule("O_2", 5.3562e-26, 213e-12, False)
CH4 = Molecule("CH_4", 2.6781e-26, 234e-12, False)
CO = Molecule("CO", 4.6867e-26, 262e-12, False)
H2 = Molecule("H_2", 3.3476e-27, 157e-12, False)

fuel = ChemReactions([O2, CH4], [CO, H2], [1, 1], [1, 2], 1e-25, 400)
print(fuel.reaction_set)
p1 = Particle(O2, [0, 0], [5, 1])
p2 = Particle(CO, [0, 0], [400, 1])
out = fuel.activation([p1, p2])
print(p1.pos)
fuel.print()
