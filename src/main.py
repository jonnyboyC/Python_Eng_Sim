from Wall import Wall
from Quadtree import Quadtree, Rectangle
from Animate import Animate
from Particle import Molecule, Particle, ChemicalReaction
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

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
    Initialize a chemical reaction based on data read from file
    :param reactants: human readable string representing the reactants
    :param products: human readable string representign the products
    :param ea: Activation energy of the reaction
    :param additional: Map that converts a chemical symbol to the molecule object
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
    Initialize a wall based on data read from file
    :param vert: Points defining the vertices of the wall
    :param rest_temp: What is the ambient temperature of the wall
    :param conduct: Thermal conductivity of the wall
    :param cp: Specific heat of the wall
    :param rho: Density of the wall
    :param d: Thickness of the wall
    :param area: Area of the wall
    :return:
    """

    vertices = re.findall("\-?[0-9]*\.[0-9]*", vert)
    for i in range(0, len(vertices)):
        vertices[i] = float(vertices[i])

    coordinates = [vertices[0:2], vertices[2:4]]
    return Wall(coordinates, rest_temp, conduct, cp, rho, d, area)


def initial_chemical(symbol, mass, radius, enthalpy, inert):
    """
    Initialize a molecule based on data read from a file
    :param symbol: Chemical symbol
    :param mass: mass of the molecule
    :param radius: hard shell radius of the molecule
    :param enthalpy: enthalpy of formation molecule
    :param inert: boolean is the molecule inert
    :return:
    """

    return Molecule(symbol, mass, radius, enthalpy, inert == 'True')

def random_start(N: int, u: float, x: float, molecules: list):
    """
    Currently used for development purposes, generates a list of particles with random, start, velocity, and
    molecule
    :param N: Number of particles to generate
    :param u: max start velocity
    :param x: max start position
    :param molecules: list of molecules that the random start can take
    :return: list of particles with random properties assigned as well as the bounded empty quaddtree
    """

    x0 = np.multiply(np.random.rand(N), 2*x)
    x1 = np.multiply(np.random.rand(N), 2*x)
    u0 = np.multiply(np.random.rand(N), 2*u)
    u1 = np.multiply(np.random.rand(N), 2*u)
    m = np.random.random_integers(0, 13, N)

    particles = []
    for i in range(0, N):
        particles.append(Particle(molecules[m[i]], [x0[i] - x, x1[i] - x], [u0[i] - u ,u1[i] - u]))

    quad = Quadtree(0, Rectangle(-1.1*x, -1.1*x, 2.2*x, 2.2*x))

    return particles, quad


def main():

    moles = r"C:\Users\John\Documents\Visual_Studio_2015\Projects\Python_Eng_Sim\test\molecules.mol"
    react = r"C:\Users\John\Documents\Visual_Studio_2015\Projects\Python_Eng_Sim\test\reactions.rec"
    walls = r"C:\Users\John\Documents\Visual_Studio_2015\Projects\Python_Eng_Sim\test\walls.wal"

    molecules = scan_file(moles, initial_chemical)
    walls = scan_file(walls, initial_walls)
    reactions = scan_file(react, initial_reaction, add_func=make_mole_map, more_input=molecules)

    N = 50
    u = 10
    x = 5e-7

    particles, quad = random_start(N, u, x, molecules)
    fig = plt.figure()

    ani = Animate(fig, walls)

    objects = []
    objects.extend(particles)
    objects.extend(walls)

    dt = 1e-10

    animate = animation.FuncAnimation(
        fig, ani, init_func=ani.init, interval=0, frames=10, blit=True,
        fargs=[particles, objects, quad, dt]
    )

    plt.show()

if __name__ == '__main__':
    main()