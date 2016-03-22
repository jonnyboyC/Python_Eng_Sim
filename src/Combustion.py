import Particle
import Wall
import numpy as np
import re
import Animate
import Quadtree
import time



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

    return Particle.ChemicalReaction(reactants, products, nu_r, nu_p, ea)


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
    return Wall.Wall(coordinates, rest_temp, conduct, cp, rho, d, area)


def initial_chemical(symbol, mass, radius, enthalpy, inert):
    """

    :param symbol:
    :param mass:
    :param radius:
    :param enthalpy:
    :param inert:
    :return:
    """

    return Particle.Molecule(symbol, mass, radius, enthalpy, inert == 'True')


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

def stub():
    pass

def collision(object0, object1):
    """
    Determines the resulting velocity after to particles elastically collide
    :param particles: list of particles
    :return: updated particles
    """

    if isinstance(object0, Wall.Wall) and isinstance(object1, Wall.Wall):
        return

    if isinstance(object0, Wall.Wall):
        inelastic_wall(wall=object0, particle=object1)

    elif isinstance(object1, Wall.Wall):
        inelastic_wall(wall=object1, particle=object0)

    else:
        inelastic_particles(object0, object1) 


def inelastic_wall(wall: Wall.Wall, particle: Particle.Particle):
    """
    return (Math.abs((l2.lat() - l1.lat())*c.lat() +  c.lng()*(l1.lng() -     
       l2.lng()) + (l1.lat() - l2.lat())*l1.lng() +
       (l1.lng() - l2.lng())*c.lat())/ Math.sqrt((l2.lat() - l1.lat())^2 +
       (l1.lng() - l2.lng())^2) <= r)
    """
    x0 = particle.x
    x1 = wall.vert[:, 1]
    unit = wall.unit

    proj_vec = (x1 - x0) - np.multiply(np.dot((x1 - x0), unit),unit)
    dis = np.linalg.norm(proj_vec) 

    if dis > particle.mole.radius:
        return

    dir_mag = np.dot(u0, proj_vec)/np.power(dis, 2)
    particle.u -= np.multiply(dir_mag, proj_vec)
    

def inelastic_particles(particle0: Particle.Particle, particle1: Particle.Particle):
    if (np.linalg.norm(particle0.x - particle1.x) > 
        particle0.mole.radius + particle1.mole.radius):
        return

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
    particle0.vel = mass_ratio0 * dir_mag * dir_vec
    particle1.vel = mass_ratio1 * dir_mag * -dir_vec

def main():

    moles = r"C:\Users\John\Documents\Visual_Studio_2015\Projects\Python_Eng_Sim\test\molecules.mol"
    react = r"C:\Users\John\Documents\Visual_Studio_2015\Projects\Python_Eng_Sim\test\reactions.rec"
    walls = r"C:\Users\John\Documents\Visual_Studio_2015\Projects\Python_Eng_Sim\test\walls.wal"

    molecules = scan_file(moles, initial_chemical)
    walls = scan_file(walls, initial_walls)
    reactions = scan_file(react, initial_reaction, add_func=make_mole_map, more_input=molecules)

    
    animation = Animate.Animate()
    animation.add_walls(walls)
    animation.expand_bounds()
    animation.show()

    N = 100
    x = np.multiply(np.random.rand(N), 0.1)
    y = np.multiply(np.random.rand(N), 0.1)
    u = np.multiply(np.random.rand(N), 20)
    v = np.multiply(np.random.rand(N), 20)
    m = np.random.random_integers(0, 13, N)

    particles = []
    quad = Quadtree.Quadtree(0, Quadtree.Rectangle(-0.06, -0.06, 0.12, 0.12))
    for i in range(0, N):
        particles.append(Particle.Particle(molecules[m[i]], [x[i]-0.05, y[i]-0.05], [u[i]-10,v[i]-10]))

    neighbors = []
    objects = []
    objects.extend(particles)
    objects.extend(walls)
    dt = 0.0001

    for i in range(1000):
        start1 = time.clock()
        for object in objects:
            quad.insert(object)
        stop1 = time.clock()

        start2 = time.clock()
        for object in objects:
            neighbors = quad.retreive(neighbors, object)
            for neighbor in neighbors:
                if object is neighbor:
                    continue
                collision(object, neighbor)
            neighbors = []
        stop2 = time.clock()

        animation.add_particles(particles)

        start3 = time.clock()
        for particle in particles:
            particle.update_x(i*dt)
        stop3 = time.clock()

        start4 = time.clock()
        quad.clear()
        stop4 = time.clock()

        print('Loop: ' + str(i) + ' insert: ' + str(stop1 - start1) + ' collision: ' + str(stop2 - start2) + ' update: ' + str(stop3 - start3) + ' clear: ' + str(stop4 - start4))

    print("you win!")

main()
