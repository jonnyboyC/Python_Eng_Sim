from Particle import Particle, Molecule, ChemicalReaction
from Wall import Wall
import numpy as np
from scipy.spatial.distance import pdist, squareform
from Quadtree import Quadtree, Rectangle


def reaction(reactions: list, particles: list):
    # Determine if collision, may result in a reaction

    # Get reactants
    """

    :type particles: object
    """
    reactants = set()
    for particle in particles:
        reactants.add(particle.symbol)

        # Check if we have an inert gas in reaction
        if particle.inert:
            collision(particles)
            return particles

    # Determine if reactants form a real reaction
    for reaction in reactions:
        if reactants >= reaction.reaction_set():
            particles = reaction.activation(particles)
            particles = collision(particles)
            return particles

    # If no reaction is present determine collision
    particles = collision(particles)
    return particles


def step_time(particles: list, objects: list, quad: Quadtree, dt: float):
    """
    TODO needs additional checks
    Step the system forward in time.
    :param particles: a list of particles currently in the system
    :param objects: a list of all particles and wall in the system
    :param quad: an empty quadtree bounding the simulation
    :param dt: timestep of the system
    """

    # Generate quadtree
    for obj in objects:
        quad.insert(obj)
            
    quad.set_leaf()

    # Detect collisions
    collision(particles, quad, dt)

    for obj in objects:
        obj.checked = False
        obj.leaf = False

    # Update positions of all particles
    for particle in particles:
        particle.update_x(dt)

    # Reset quadtree
    quad.clear()


def collision(particles: list, quad: Quadtree, dt: float):
    """
    Determine if particles collide and update velocities according to colission type
    :param particles: list of particles in the simulation
    :param quad: Quadtree containing particles
    :param dt: timestep in seconds
    """
    for particle in particles:
        if particle.checked:
            continue

        neighbors = []

        # Retrieve objects that may intersect this particle
        neighbors = quad.retreive(neighbors, particle)

        # separate particles from walls
        p_neighbors = [neighbor for neighbor in neighbors 
                        if isinstance(neighbor, Particle) 
                        and neighbor.checked is False]
        w_neighbors = [neighbor for neighbor in neighbors
                        if isinstance(neighbor, Wall)]

        # Determine which objects at a given level of the quadtree itersect
        positions = np.array([neighbor.x for neighbor in p_neighbors])
        distances = squareform(pdist(positions))
        radius_sums = radius_sum(p_neighbors)

        idx1, idx2 = np.where(distances - radius_sums < 0)
        unique = (idx1 < idx2)

        idx1 = idx1[unique]
        idx2 = idx2[unique]

        # If two particles are shown to collide update particles
        for id1, id2 in zip(idx1, idx2):
            inelastic_particles(p_neighbors[id1], p_neighbors[id2], dt)
            particles[id1].checked = True
            particles[id2].checked = True

        # Check particles from this retreive to walls retreived
        for wall in w_neighbors:
            for particle in p_neighbors:
                if inelastic_wall(wall, particle, dt):
                    particle.checked = True

        # If particle is a leaf node it has been checked againsted all canidates
        for neighbor in p_neighbors:
            if neighbor.leaf:
                neighbor.checked = True


def inelastic_wall(wall: Wall, particle: Particle, dt: float):
    """
    TODO extending for moveable walls
    Determines the resulting velocity of a particle bouncing off a wall, with mid-frame position correction
    :param wall: the wall object that may intersect the particle
    :param particle: the particle in the vicinity of the wall
    :param dt: length of the timestep
    :return: Boolean value of true if a collision occured
    """
    x0 = particle.x
    u0 = particle.u
    x1 = wall.vert[0, :]
    p = x0 - x1
    r = particle.mole.radius

    unit = wall.unit
    norm = wall.norm
    dis = np.dot(p, norm)

    # Determine if a collision actually occured
    if dis > r:
        return False


    #if norm_dir > 0:
    #    return False

    # Decompose particles velocity normal and tangential to wall
    norm_dir = np.dot(u0, norm)
    nu = np.multiply(norm_dir, norm)
    uu = np.multiply(np.dot(u0, unit), unit)

    # Update velocity and position
    x0 -= u0 * dt
    dt0 = (np.abs(np.dot(x0 - x1, norm)) - r) / np.linalg.norm(nu)
    x0 += u0 * dt0
    particle.u = uu - nu
    particle.x = x0 + particle.u*(dt - dt0)

    return True
    

def inelastic_particles(particle0: Particle, particle1: Particle, dt: float):
    """
    TODO extending for moveable walls
    Determines the resulting velocity of a set of colliding particles with mid-fram correction
    :param particle0: First particle in the collision
    :param particle1: Second particle in the collision
    """

    # Determine mass ratios
    mass_ratio0 = 2 * particle1.mole.mass / (particle0.mole.mass + particle1.mole.mass)
    mass_ratio1 = 2 * particle0.mole.mass / (particle0.mole.mass + particle1.mole.mass)

    # generate relative terms
    X = particle0.x - particle1.x
    U = particle1.u - particle0.u
    R = particle0.mole.radius + particle1.mole.radius

    a = np.dot(U, U)
    b = 2*np.dot(X, U)
    c = np.dot(X, X) - R*R

    # Determine time of collision
    dt0 = (-b + np.sqrt(b*b - 4*a*c))/(2*a)

    x0 = particle0.x - particle0.u*dt0
    x1 = particle1.x - particle1.u*dt0

    u0 = particle0.u
    u1 = particle1.u

    # Direction magnitude
    dir_mag = np.dot(u0 - u1, x0 - x1) / (np.dot(x0 - x1, x0 - x1))
    dir_vec = (x0 - x1)

    # Set new velocities
    particle0.u -= mass_ratio0 * dir_mag * dir_vec
    particle1.u -= mass_ratio1 * dir_mag * -dir_vec

    # Correct positions
    particle0.x = x0 + particle0.u * dt0
    particle1.x = x1 + particle1.u * dt0


def radius_sum(particles: list):
    """
    Determines each pairwise sum of elements in particles for intersection
    :param particles: the particle currently being checked
    :return np.array: size [len(particles), len(particles)] of pairwise sums
    """
    radius = np.array([particle.mole.radius for particle in particles])
    sums = radius[:, None] + radius[None, :]

    return sums
