"""Module of axes getter functions."""

import pathlib
from itertools import combinations

import numpy as np
import pywindow as pw
import rdkit.Chem.AllChem as rdkit  # noqa: N813
import stk
from scipy.spatial.distance import pdist
from sklearn.metrics.pairwise import cosine_similarity

from .utils import distance, normalize_vector


def by_pywindow(filename: pathlib.Path | str) -> np.ndarray:
    """Get axes using pyWindow windows."""
    # This doesnt work, not sure I understand pywindow
    molsys = pw.MolecularSystem.load_file(filename)
    mol = molsys.system_to_molecule()
    windows = mol.calculate_windows()
    com = mol.calculate_centre_of_mass()
    return windows - com


def by_smarts(molecule: stk.Molecule, smarts_string: str) -> np.ndarray:
    """Get axes using smarts strings of a subgroup."""
    # Not finished?     # return "Result of BySmarts"
    raise NotImplementedError
    rdkit_mol = molecule.to_rdkit_mol()
    rdkit.SanitizeMol(rdkit_mol)
    centroid_smiles = [
        molecule.get_centroid(atom_ids=atom_ids)
        for atom_ids in rdkit_mol.GetSubstructMatches(
            query=rdkit.MolFromSmarts(smarts_string)
        )
    ]
    centroid_smiles = np.asarray(centroid_smiles)
    centroid_mol = molecule.get_centroid()
    distances = [
        np.linalg.norm(smile - centroid_mol) for smile in centroid_smiles
    ]
    vectors = np.array(
        [
            (smile - centroid_mol) / np.linalg.norm(smile - centroid_mol)
            for smile in centroid_smiles
        ]
    )
    return vectors, np.mean(distances)


def by_smiles(
    molecule: stk.Molecule,
    smiles_string: str,
) -> tuple[np.ndarray, float]:
    """Get axes using smiles strings of a subgroup."""
    rdkit_mol = molecule.to_rdkit_mol()
    rdkit.SanitizeMol(rdkit_mol)
    centroid_smiles = [
        molecule.get_centroid(atom_ids=atom_ids)
        for atom_ids in rdkit_mol.GetSubstructMatches(
            query=rdkit.MolFromSmiles(smiles_string)
        )
    ]
    centroid_smiles = np.asarray(centroid_smiles)
    centroid_mol = molecule.get_centroid()
    distances = [
        np.linalg.norm(smile - centroid_mol) for smile in centroid_smiles
    ]
    vectors = np.array(
        [
            (smile - centroid_mol) / np.linalg.norm(smile - centroid_mol)
            for smile in centroid_smiles
        ]
    )
    return vectors, np.mean(distances)


def by_midpoint(
    vectors: np.ndarray,
    vertice_size: float,
    no_vectors_define_facet: int,
    tolerance: float = 0.1,
) -> np.ndarray:
    """Get axes using the midpoint of given vectors."""
    if isinstance(vectors, list) and isinstance(vectors[0], np.ndarray):
        vectors = vectors[0]

    # Convert numpy array of vectors into a list of numpy arrays
    vectors_list = [vectors[i] for i in range(vectors.shape[0])]
    all_distances = [
        distance(v1, v2) for v1, v2 in combinations(vectors_list, 2)
    ]
    min_distance = min(filter(lambda d: d > 0, all_distances))

    midpoints = []  # List to store all midpoints that meet the condition
    midpoint_size = []  # List to store all midpoints that meet the condition

    # Check each combination of n vectors
    for combination in combinations(vectors_list, no_vectors_define_facet):
        distances = [
            distance(
                combination[i],
                combination[(i + 1) % no_vectors_define_facet],
            )
            for i in range(no_vectors_define_facet)
        ]
        if max(distances) - min(distances) <= tolerance * min_distance:
            midpoint = np.mean(combination, axis=0)
            midpoint_size.append(np.linalg.norm(midpoint))
            midpoint = normalize_vector(midpoint)
            midpoints.append(midpoint)  # Add the computed midpoint to the list

    return np.array(midpoints), np.mean(midpoint_size) * vertice_size


def by_midpoint_mindistance(
    vectors: np.ndarray,
    vertice_size: float,
    no_vectors_define_facet: int,
    tolerance: float = 0.1,
) -> np.ndarray:
    """Get axes using the midpoint of given vectors with a min-distance."""
    if isinstance(vectors, list) and isinstance(vectors[0], np.ndarray):
        vectors = vectors[0]

    # Convert numpy array of vectors into a list of numpy arrays
    vectors_list = [vectors[i] for i in range(vectors.shape[0])]
    all_distances = [
        distance(v1, v2) for v1, v2 in combinations(vectors_list, 2)
    ]
    min_distance = min(filter(lambda d: d > 0, all_distances))

    midpoints = []
    midpoint_size = []

    # Check each combination of n vectors
    for combination in combinations(vectors_list, no_vectors_define_facet):
        distances = [
            distance(
                combination[i],
                combination[(i + 1) % no_vectors_define_facet],
            )
            for i in range(no_vectors_define_facet)
        ]
        if np.isclose(
            max(distances), min(distances), atol=tolerance, rtol=0
        ) and np.isclose(max(distances), min_distance, atol=tolerance, rtol=0):
            midpoint = np.mean(combination, axis=0)
            midpoint_size.append(np.linalg.norm(midpoint))
            midpoint = normalize_vector(midpoint)
            midpoints.append(midpoint)

    return np.array(midpoints), np.mean(midpoint_size) * vertice_size


def by_midpoint_square(
    vectors: np.ndarray,
    vertice_size: float,
    no_vectors_define_facet: int,
    tolerance: float = 0.1,
) -> np.ndarray:
    """Get axes using the midpoint of rectangular vectors."""
    if isinstance(vectors, list) and isinstance(vectors[0], np.ndarray):
        vectors = vectors[0]

    # Convert numpy array of vectors into a list of numpy arrays
    vectors_list = [vectors[i] for i in range(vectors.shape[0])]
    all_distances = [
        distance(v1, v2) for v1, v2 in combinations(vectors_list, 2)
    ]
    min_distance = min(filter(lambda d: d > 0, all_distances))

    midpoints = []  # List to store all midpoints that meet the condition
    midpoint_size = []  # List to store all midpoints that meet the condition

    # Check each combination of n vectors
    for combination in combinations(vectors_list, no_vectors_define_facet):
        distances = sorted(pdist(combination))
        min_found_distance = min(distances)
        max_found_distance = max(distances)
        if (
            np.isclose(
                min_found_distance, min_distance, atol=tolerance, rtol=0
            )
            and np.allclose(
                [max_found_distance], distances[4:], atol=tolerance, rtol=0
            )
            and np.allclose(
                [min_found_distance], distances[:4], atol=tolerance, rtol=0
            )
        ):
            midpoint = np.mean(combination, axis=0)
            midpoint_size.append(np.linalg.norm(midpoint))
            midpoint = normalize_vector(midpoint)
            midpoints.append(midpoint)  # Add the computed midpoint to the list

    return np.array(midpoints), np.mean(midpoint_size) * vertice_size


def remove_common(
    arr1: np.ndarray,
    arr2: np.ndarray,
    tolerance: float = 0.1,
) -> np.ndarray:
    """Remove common vectors."""
    filtered_arr1 = []
    for vec1 in arr1:
        diff = True
        for _, vec2 in enumerate(arr2):
            similarity = cosine_similarity([vec1], [vec2])[0][0]
            if similarity > (1 - tolerance):
                diff = False
                break
        if diff:
            filtered_arr1.append(vec1)

    return np.array(filtered_arr1)
