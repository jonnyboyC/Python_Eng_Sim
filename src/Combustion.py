import Particle

def valid_reaction(reactions: list, particles: list):
    # Determine if collision, may result in a reaction

    # Get reactants
    reactants = set()
    for particle in particles:
        reactants.add(particle.symbol)

        # Check if we have an inert gas in reaction
        if particle.inert:
            particles = collision(particles)
            return particles

    # Determine if reactants form a real reaction
    for reaction in reactions:
        if reactants == reaction.reaction_set():
            particles = reaction.activation(particles):
            return particles

    # If no reaction is present determine collision
    particles = collision(particles)
    return particles


def conditions(particles: list):
    # Check if conditions are right for reaction

    """

    :type particles: list
    """
    if True:
        particles[0] = collision(particles[0])
        return particles


def collision(particles: list):

    # Determine mass ratios
    mass_ratio0 = 2*particles[1].mass/(particles[0].mass + particles[1].mass)
    mass_ratio1 = 2*particles[0].mass/(particles[0].mass + particles[1].mass)

    # Get position and velocity vectors
    v0 = particles[0].vel
    v1 = particles[1].vel
    x0 = particles[0].pos
    x1 = particles[1].pos

    # Direction magnitude
    dir_mag = np.dot(v0-v1, x0-x1)/(np.power(np.linalg.norm(x0-x1), 2))
    dir_vec = (x0 - x1)

    # Set new velocities
    particles[0].vel = mass_ratio0*dir_mag*dir_vec
    particles[1].vel = mass_ratio1*dir_mag*-dir_vec

    return particles