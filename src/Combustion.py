from Particle import Particle, Molecule, ChemicalReaction
from Wall import Wall
import numpy as np
from scipy.spatial.distance import pdist, squareform
import re
import Animate 
from  Quadtree import Quadtree, Rectangle
import time
import matplotlib



def scan_file(path: str, initial_func, add_func=None, more_input=None):
    """
    Generic function for reading text files to initialize simulation objects
    :param path: Path of the file containing object information
    :param initial_func: function pointer taking parameters gathered and return an initialized object
    :param add_func: function pointer taking input list more_input returning an output
    :param more_input: list of inputs used by function add_func
    :return: list containing the initialized objects
    """

    parameters = []
    output = []
    file = open(path)

    if add_func is not None:
        more_input = add_func(more_input)

    for line in file:
        if line[0] == '#':
            continue

        if line[0] == '\n':
            if len(parameters) == 0:
                break

            if more_input is None:
                output.append(
                    initial_func(*parameters)
                )
            else:
                output.append(
                    initial_func(*parameters, additional=more_input)
                )

            parameters = []
            continue

        parameters.append(line.replace('\n', ''))

    return output


def make_mole_map(molecules: list):
    """
    Take a list of molecule objects and return a map relating chemical symbol to object
    :param molecules: list of molecules
    :return: symbol map to object
    """

    mole_map = {}
    for molecule in molecules:
        mole_map[molecule.symbol] = molecule

    return mole_map

def initial_reaction(reactants, products, ea, additional):
    """

    :param reactants:
    :param products:
    :param ea:
    :param additional:
    :return:
    """

    reactants = re.findall('[0-9][0-Z]+', reactants)
    products = re.findall('[0-9][0-Z]+', products)

    nu_r = []
    nu_p = []

    for react in reactants:
        nu_r.append(int(react[0]))

    for prod in products:
        nu_p.append(int(prod[0]))

    reactants[:] = [additional[a[1:]] for a in reactants]
    products[:] = [additional[a[1:]] for a in products]

    return ChemicalReaction(reactants, products, nu_r, nu_p, ea)


def initial_walls(vert, rest_temp, conduct, cp, rho, d, area):
    """

    :param vert:
    :param rest_temp:
    :param conduct:
    :param cp:
    :param rho:
    :param d:
    :param area:
    :return:
    """

    vertices = re.findall("\-?[0-9]*\.[0-9]*", vert)
    for i in range(0, len(vertices)):
        vertices[i] = float(vertices[i])

    coordinates = [vertices[0:2], vertices[2:4]]
    return Wall(coordinates, rest_temp, conduct, cp, rho, d, area)


def initial_chemical(symbol, mass, radius, enthalpy, inert):
    """

    :param symbol:
    :param mass:
    :param radius:
    :param enthalpy:
    :param inert:
    :return:
    """

    return Molecule(symbol, mass, radius, enthalpy, inert == 'True')


def reaction(reactions: list, particles: list):
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
        if reactants >= reaction.reaction_set():
            particles = reaction.activation(particles)
            particles = collision(particles)
            return particles

    # If no reaction is present determine collision
    particles = collision(particles)
    return particles


def step_time(i: int, particles: list, objects: list, quad: Quadtree, dt: float): 
    #start1 = time.clock()
    for object in objects:
            quad.insert(object)
            
    quad.set_leaf()
    #stop1 = time.clock()
  
    #start2 = time.clock()
    collision(particles, quad)

    for object in objects:
        object.checked = False
        object.leaf = False
    #stop2 = time.clock()

    #start3 = time.clock()
    for particle in particles:
        particle.update_x(dt)
    #stop3 = time.clock()

    #start4 = time.clock()
    quad.clear()
    #stop4 = time.clock()

    #print('Loop: ' + str(i) + ' insert: ' + str(stop1 - start1) + ' collision: ' + str(stop2 - start2) + ' update: ' + str(stop3 - start3) + ' clear: ' + str(stop4 - start4))

def collision(particles: list, quad: Quadtree):
    """
    Determine if particles collide and update velocities according to colission type
    :param particles: list of particles in the simulation
    :param quad: Quadtree containing particles
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
            inelastic_particles(particles[id1], particles[id2])
            particles[id1].checked = True
            particles[id2].checked = True

        # Check particles from this retreive to walls retreived
        for wall in w_neighbors:
            for particle in p_neighbors:
                if inelastic_wall(wall, particle):
                    particle.checked = True

        # If particle is a leaf node it has been checked againsted all canidates
        for neighbor in p_neighbors:
            if neighbor.leaf:
                neighbor.checked = True


def inelastic_wall(wall: Wall, particle: Particle):
    """
    TODO extending for moveable walls
    Determines the resulting velocity of a particle bouncing off a wall
    :param wall: the wall object that may intersect the particle
    :param particle: the particle in the vicinity of the wall
    :return: Boolean value of true if a collision occured
    """
    x0 = particle.x
    u0 = particle.u
    x1 = wall.vert[:, 1]
    p = x0 - x1
    unit = wall.unit
    norm = wall.norm

    dis = np.abs(np.dot(p, norm))

    if dis > particle.mole.radius:
        return False

    norm_dir = np.dot(u0, norm)
    if norm_dir > 0:
        return

    nu = np.multiply(norm_dir, norm)
    uu = np.multiply(np.dot(u0, unit), unit)
    particle.u = uu - nu

    return True
    

def inelastic_particles(particle0: Particle, particle1: Particle):
    """
    TODO extending for moveable walls
    Determines the resulting velocity of a particle bouncing off a wall
    :param wall: the wall object that may intersect the particle
    :param particle: the particle in the vicinity of the wall
    :return: Boolean value of true if a collision occured
    """

    # Determine mass ratios
    mass_ratio0 = 2 * particle1.mole.mass / (particle0.mole.mass + particle1.mole.mass)
    mass_ratio1 = 2 * particle0.mole.mass / (particle0.mole.mass + particle1.mole.mass)

    # Get position and velocity vectors
    u0 = particle0.u
    u1 = particle1.u
    x0 = particle0.x
    x1 = particle1.x

    # Direction magnitude
    dir_mag = np.dot(u0 - u1, x0 - x1) / (np.power(np.linalg.norm(x0 - x1), 2))
    dir_vec = (x0 - x1)

    # Set new velocities
    particle0.u += mass_ratio0 * dir_mag * dir_vec
    particle1.u += mass_ratio1 * dir_mag * -dir_vec


def radius_sum(particles: list):
    """
    Determines each pairwise sum of elements in particles for intersection
    :param particles: the particle currently being checked
    :return np.array: size [len(particles), len(particles)] of pairwise sums
    """
    radius = np.array([particle.mole.radius for particle in particles])
    sums = radius[:, None] + radius[None, :]

    return sums


def random_start(N: int, u: float, x: float, molecules: list):
    """

    :param N:
    :param u:
    :param x:
    :param molecules:
    :return:
    """

    x0 = np.multiply(np.random.rand(N), 2*x)
    x1 = np.multiply(np.random.rand(N), 2*x)
    u0 = np.multiply(np.random.rand(N), 2*u)
    u1 = np.multiply(np.random.rand(N), 2*u)
    m = np.random.random_integers(0, 13, N)

    particles = []
    for i in range(0, N):
        particles.append(Particle(molecules[m[i]], [x0[i] - x, x1[i] - x], [u0[i] - u ,u1[i] - u]))

    quad = Quadtree(0, Rectangle(-1.2*x, 1.2*x, 2.4*x, 2.4*x))

    return particles, quad


def main():

    moles = r"C:\Users\John\Documents\Visual_Studio_2015\Projects\Python_Eng_Sim\test\molecules.mol"
    react = r"C:\Users\John\Documents\Visual_Studio_2015\Projects\Python_Eng_Sim\test\reactions.rec"
    walls = r"C:\Users\John\Documents\Visual_Studio_2015\Projects\Python_Eng_Sim\test\walls.wal"

    molecules = scan_file(moles, initial_chemical)
    walls = scan_file(walls, initial_walls)
    reactions = scan_file(react, initial_reaction, add_func=make_mole_map, more_input=molecules)

    particles, quad = random_start(1000, 10,  5e-7, molecules)
    # fig = Animate.make_fig()

    # ani = matplotlib.animation.FuncAnimation(fig, Animate.animate, np.arange(0, 200), 

    objects = []
    objects.extend(particles)
    objects.extend(walls)
    dt = 1e-10

    for i in range(1000):
        step_time(i, particles, objects, quad, dt)
    
    # matplotlib.animation.FuncAnimation(
    print("you win!")


if __name__ == '__main__':
    main()
