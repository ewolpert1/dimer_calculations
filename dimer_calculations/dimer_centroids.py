"""Dimer centroids module."""

import logging
import pathlib
import re

import networkx as nx
import numpy as np
import stk
import stko


def parse_atom_and_bond_data(data: list[str]) -> tuple:
    """Parse .mol file data and extract atom coordinates and bond information.

    Arguments:
    ---------
        data (list of str):
            Lines from the .mol file.

    Returns:
    -------
        tuple:
            List of (x, y, z) coordinates for each atom, List of
            (atom1, atom2) bonds.

    """
    atom_data = []
    bond_data = []
    atom_section, bond_section = False, False

    for line in data:
        if "M  V30 BEGIN ATOM" in line:
            atom_section = True
            continue
        if "M  V30 END ATOM" in line:
            atom_section = False
            continue
        if "M  V30 BEGIN BOND" in line:
            bond_section = True
            continue
        if "M  V30 END BOND" in line:
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

 
def calculate_centroid(atoms: list[tuple]) -> tuple:
    """Calculate the centroid of a group of atoms.

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


def separate_molecules(atom_data: list, bond_data: list) -> list:
    """Separate the molecules on bond data and return the atoms in each.

    Args:
    ----
    atom_data (list of tuples): (index, (x, y, z)) coordinates for each atom.
    bond_data (list of tuples): (atom1, atom2) bonds.

    Returns:
    -------
    list of lists: Each sublist contains the (x, y, z) coordinates of atoms in
    one molecule.

    """
    graph = nx.Graph()
    graph.add_edges_from(bond_data)

    # Creating a mapping from atom index to coordinates
    atom_dict = dict(atom_data)

    molecules = []
    for component in nx.connected_components(graph):
        molecule = []
        for i in component:
            if atom_dict.get(i) is None:
                logging.info(f"Key {i} not found in atom_dict")
            else:
                molecule.append(atom_dict[i])
        molecules.append(molecule)

    return molecules


def return_centroids(mol_file: str | pathlib.Path) -> list[np.ndarray]:
    """Get centroids of mol file."""
    # Read the .mol file
    with pathlib.Path(mol_file).open() as f:
        mol_data = f.readlines()

    # Parsing atom and bond data
    atom_data, bond_data = parse_atom_and_bond_data(mol_data)

    # Separating the molecules
    molecules = separate_molecules(atom_data, bond_data)

    # Calculating centroids for each molecule
    return [calculate_centroid(molecule) for molecule in molecules]


def get_dimer_centroids(
    dimer_molecule: stk.Molecule,
) -> tuple[np.ndarray, np.ndarray]:
    """Get centroids of dimers."""
    graph = stko.Network.init_from_molecule(dimer_molecule)

    # A series of graphs still connected.
    connected_graphs = graph.get_connected_components()

    centroids = []
    for cg in connected_graphs:
        # Get atoms from nodes.
        atoms = list(cg)
        atom_ids = tuple(i.get_id() for i in atoms)

        position_matrix = np.array(
            tuple(
                i
                for i in dimer_molecule.get_atomic_positions(atom_ids=atom_ids)
            )
        )

        centroid = np.divide(
            position_matrix.sum(axis=0),
            len(atom_ids),
        )
        centroids.append(centroid)
    return centroids
