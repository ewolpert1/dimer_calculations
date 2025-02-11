import re

import networkx as nx


def parse_atom_and_bond_data(data):
    """Parses the .mol file data and extracts atom coordinates and bond information.

    Args:
    ----
    data (list of str): Lines from the .mol file.

    Returns:
    -------
    tuple: List of (x, y, z) coordinates for each atom, List of (atom1, atom2) bonds.

    """
    atom_data = []
    bond_data = []
    atom_section, bond_section = False, False

    for line in data:
        if "M  V30 BEGIN ATOM" in line:
            atom_section = True
            continue
        elif "M  V30 END ATOM" in line:
            atom_section = False
            continue
        elif "M  V30 BEGIN BOND" in line:
            bond_section = True
            continue
        elif "M  V30 END BOND" in line:
            break

        if atom_section:
            # Extracting the X, Y, Z coordinates from each atom line
            match = re.search(
                r"M  V30 (\d+) \w+ (-?\d+\.\d+) (-?\d+\.\d+) (-?\d+\.\d+)",
                line,
            )
            if match:
                atom_index, x, y, z = match.groups()
                atom_data.append(
                    (int(atom_index), (float(x), float(y), float(z)))
                )
        elif bond_section:
            # Extracting the bond information
            match = re.search(r"M  V30 \d+ \d+ (\d+) (\d+)", line)
            if match:
                atom1, atom2 = map(int, match.groups())
                bond_data.append((atom1, atom2))

    return atom_data, bond_data


def calculate_centroid(atoms):
    """Calculates the centroid of a group of atoms.

    Args:
    ----
    atoms (list of tuples): (x, y, z) coordinates for each atom.

    Returns:
    -------
    tuple: (x, y, z) coordinates of the centroid.

    """
    x_sum, y_sum, z_sum = 0, 0, 0
    for x, y, z in atoms:
        x_sum += x
        y_sum += y
        z_sum += z
    count = len(atoms)
    return (x_sum / count, y_sum / count, z_sum / count)


def separate_molecules(atom_data, bond_data):
    """Separates the molecules based on bond data and returns the atoms in each molecule.

    Args:
    ----
    atom_data (list of tuples): (index, (x, y, z)) coordinates for each atom.
    bond_data (list of tuples): (atom1, atom2) bonds.

    Returns:
    -------
    list of lists: Each sublist contains the (x, y, z) coordinates of atoms in one molecule.

    """
    G = nx.Graph()
    G.add_edges_from(bond_data)

    # Creating a mapping from atom index to coordinates
    atom_dict = dict(atom_data)

    molecules = []
    for component in nx.connected_components(G):
        # print(component, len(component),len(atom_dict),atom_dict[2])
        # molecule = [atom_dict[i] for i in component]
        molecule = []
        for i in component:
            if atom_dict.get(i) is None:
                print(f"Key {i} not found in atom_dict")
            else:
                molecule.append(atom_dict[i])
        molecules.append(molecule)

    return molecules


def return_centroids(mol_file):
    # Read the .mol file
    with open(mol_file) as file:
        mol_data = file.readlines()

    # Parsing atom and bond data
    atom_data, bond_data = parse_atom_and_bond_data(mol_data)

    # Separating the molecules
    molecules = separate_molecules(atom_data, bond_data)

    # Calculating centroids for each molecule
    centroids = [calculate_centroid(molecule) for molecule in molecules]

    return centroids
