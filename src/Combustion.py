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


def collision(particles: list):
    """
    Determines the resulting velocity after to particles elastically collide
    :param particles: list of particles
    :return: updated particles
    """

    # Determine mass ratios
    mass_ratio0 = 2 * particles[1].mass / (particles[0].mass + particles[1].mass)
    mass_ratio1 = 2 * particles[0].mass / (particles[0].mass + particles[1].mass)

    # Get position and velocity vectors
    v0 = particles[0].vel
    v1 = particles[1].vel
    x0 = particles[0].pos
    x1 = particles[1].pos

    # Direction magnitude
    dir_mag = np.dot(v0 - v1, x0 - x1) / (np.power(np.linalg.norm(x0 - x1), 2))
    dir_vec = (x0 - x1)

    # Set new velocities
    particles[0].vel = mass_ratio0 * dir_mag * dir_vec
    particles[1].vel = mass_ratio1 * dir_mag * -dir_vec

    return particles


# path = input("raise your dong: ")
moles = r"C:\Users\John\PycharmProjects\Python_Eng_Sim\test\molecules.mol"
react = r"C:\Users\John\PycharmProjects\Python_Eng_Sim\test\reactions.rec"
walls = r"C:\Users\John\PycharmProjects\Python_Eng_Sim\test\walls.wal"

molecules = scan_file(moles, initial_chemical)
walls = scan_file(walls, initial_walls)
reactions = scan_file(react, initial_reaction, add_func=make_mole_map, more_input=molecules)

animation = Animate.Animate()
animation.show()
animation.add_walls(walls)
animation.expand_bounds()

N = 500
x = np.random.rand(N)
y = np.random.rand(N)
u = np.multiply(np.random.rand(N), 10)
v = np.multiply(np.random.rand(N), 10)
m = np.random.random_integers(0, 13, N)

particles = []
quad = Quadtree.Quadtree(0, Quadtree.Rectangle(-0.5, -0.5, 1, 1))
for i in range(0, N):
    particles.append(Particle.Particle(molecules[m[i]], [x[i]-0.5, y[i]-0.5], [u[i],v[i]]))

start = time.clock()
for i in range(0, 1000):
    for particle in particles:
        quad.insert(particle.aabb())
    quad.clear()
end = time.clock()
print(str((end-start)/1000))
animation.add_particles(particles)

print("you win!")
