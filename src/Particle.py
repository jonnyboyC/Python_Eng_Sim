import numpy as np


class ValidReactions:
    """
    derp
    """

    def __init__(self, file: str):
        file = open(file, 'r')
        pass
        t = []
        self.reactions = set(frozenset(i) for i in t)


class Particle:
    """
    derp
    """

    def __init__(self, symbol: str, mass: float, vel, pos, temp: float, inert: bool):
        self.symbol = symbol
        self.mass = float(mass)
        self.vel = np.array(vel)
        self.pos = np.array(pos)
        self.temp = float(temp)
        self.inert = bool(inert)

    def update_velocity(self, dt):
        # Update particles position
        self.pos += dt*self.vel

    def radiate_heat(self, dt):
        pass


def valid_reaction(valid_list: ValidReactions, particles: list):
    # Determine if collision, may result in a reaction

    # Get reactants
    reactants = {}
    for particle in particles:
        reactants.add(particle.symbol)

    # Determine if reactants are a real reaction
    for valid in valid_list.reactions:
        if reactants == valid:
            particles = valid_list.conditions(particles, valid_list)
            return particles

    # If no reaction is present determine collision
    particles = collision(particles)
    return particles


def conditions(particles: list):
    # Check if conditions are right for reaction
    if True:
        particles = collision(particles)
        return particles



def collision(particles: list):

    # Get relevant instance variables
    v1 = particles[1].vel
    v2 = particles[1].vel
    m1 = particles[1].mass
    m2 = particles[2].mass
    x1 = particles[1].pos
    x2 = particles[2].pos

    # Set new velocities assuming inelastic collisions
    particles[1].vel = \
        ((v1 - v2)*2*m2/(m1 + m2)*(x1-x2)*
        np.dot(v1-v2, x1-x2)/(np.power(np.power(np.linalg.norm(x1-x2)),2)))
    particles[2].vel = \
        ((v2 - v1)*2*m1/(m1 + m2)*(x1-x2)*
        np.dot(v2-v1, x2-x1)/(np.power(np.power(np.linalg.norm(x1-x2)),2)))

    return particles






