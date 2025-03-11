"""Module to check catenation."""

import numpy as np
import pore_mapper as pm
import stk
import stko
from scipy.spatial.distance import cdist


def two_pores(
    stk_molecule: stk.Molecule,
    centroid1: np.ndarray,
    centroid2: np.ndarray,
) -> list[float]:
    """Run poremapper on dimer."""
    # Read in host from xyz file.
    host = pm.Host(
        atoms=(
            pm.Atom(id=i.get_id(), element_string=i.__class__.__name__)
            for i in stk_molecule.get_atoms()
        ),
        position_matrix=stk_molecule.get_position_matrix(),
    )

    # Number of centroids can be any number:
    pore_distances = []
    centroids = (centroid1, centroid2)

    for centroid in centroids:
        # Define calculator object.
        calculator = pm.Inflater(bead_sigma=1.2, centroid=centroid)

        # Run calculator on host object, analysing output.
        final_result = calculator.get_inflated_blob(host=host)
        pore = final_result.pore
        pore_distances.append(pore.get_mean_distance_to_com())

    return pore_distances


def one_pore(stk_molecule: stk.Molecule) -> float:
    """Run poremapper on one cage."""
    # Read in host from xyz file.
    host = pm.Host(
        atoms=(
            pm.Atom(id=i.get_id(), element_string=i.__class__.__name__)
            for i in stk_molecule.get_atoms()
        ),
        position_matrix=stk_molecule.get_position_matrix(),
    )

    # Number of centroids can be any number:
    pore_distances = []
    calculator = pm.Inflater(bead_sigma=1.2, centroid=host.get_centroid())

    final_result = calculator.get_inflated_blob(host=host)
    pore = final_result.pore
    pore_distances.append(pore.get_mean_distance_to_com())

    return pore_distances


def two_com_distances(stk_molecule: stk.Molecule) -> list[float]:
    """Get com distances on dimer."""
    graph = stko.Network.init_from_molecule(stk_molecule)

    # A series of graphs still connected.
    connected_graphs = graph.get_connected_components()

    pore_distances = []
    for cg in connected_graphs:
        # Get atoms from nodes.
        atoms = list(cg)
        atom_ids = tuple(i.get_id() for i in atoms)

        pair_dists = cdist(
            stk_molecule.get_position_matrix(),
            stk_molecule.get_centroid(atom_ids).reshape(1, 3),
        )
        pore_distances.append(np.min(pair_dists.flatten()))

    return pore_distances


def one_com_distance(stk_molecule: stk.Molecule) -> float:
    """Get com distance on one cage."""
    return stko.molecule_analysis.GeometryAnalyser().get_min_centroid_distance(
        stk_molecule
    )


def get_centroids(stk_molecule: stk.Molecule) -> list[np.ndarray]:
    """Get centroids of molecule."""
    graph = stko.Network.init_from_molecule(stk_molecule)

    # A series of graphs still connected.
    connected_graphs = graph.get_connected_components()

    centroids = []
    for cg in connected_graphs:
        # Get atoms from nodes.
        atoms = list(cg)
        atom_ids = tuple(i.get_id() for i in atoms)
        centroids.append(stk_molecule.get_centroid(atom_ids))

    return centroids


def check_catenane(
    one_pore: float | list[float], 
    two_pore: list[float],
    cutoff: float = 0.2,
) -> bool:
    """Check if two cages are catenated."""
    if isinstance(one_pore, list):
        one_pore = one_pore[0] 
    return bool(
        two_pore[0] < (one_pore - cutoff) and two_pore[1] < (one_pore - cutoff)
    )
