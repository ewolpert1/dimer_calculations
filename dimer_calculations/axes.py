"""Module of axes getter functions."""

from itertools import combinations

import numpy as np
import rdkit.Chem.AllChem as rdkit  # noqa: N813
import stk
from scipy.spatial.distance import pdist
from sklearn.metrics.pairwise import cosine_similarity

from .utils import distance, normalize_vector


def by_smarts(molecule: stk.Molecule, smarts_string: str) -> np.ndarray:
    """Get axes using smarts strings of a subgroup."""
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


def by_mol_file(self, mol_file, removed_structure="[NX3H2]"):
    m = rdkit.MolFromMolFile(mol_file)
    diamine_mol_block = open(mol_file).read()
    diamine_mol = rdkit.MolFromMolBlock(diamine_mol_block)
    amino_smarts = (
        removed_structure  # SMARTS pattern for primary amine nitrogen
    )
    amino_query = rdkit.MolFromSmarts(amino_smarts)
    matches = diamine_mol.GetSubstructMatches(amino_query)
    rw_mol = rdkit.RWMol(diamine_mol)
    atoms_to_remove = set()
    for match in matches:
        n_idx = match[0]
        # Remove attached hydrogens
        for neighbor in rw_mol.GetAtomWithIdx(n_idx).GetNeighbors():
            if neighbor.GetAtomicNum() == 1:
                atoms_to_remove.add(neighbor.GetIdx())
        # Mark the nitrogen atom for removal
        atoms_to_remove.add(n_idx)
    # Remove atoms in reverse order to avoid reindexing issues
    for atom_idx in sorted(atoms_to_remove, reverse=True):
        rw_mol.RemoveAtom(atom_idx)
    diamine_core = rw_mol.GetMol()
    diamine_core.UpdatePropertyCache()
    # Generate the substructure query from the diamine core
    substructure_query = diamine_core
    # Perform substructure matching on the cage molecule
    matches = cage_mol.GetSubstructMatches(substructure_query)
    if not matches:
        print(f"No matches found for diamine core in entry {idx}.")
    # Calculate centroids of the matched substructures
    match_centroids = []
    for match in matches:
        atom_positions = np.array([conf.GetAtomPosition(idx) for idx in match])
        centroid = get_centroid(atom_positions)
        match_centroids.append(centroid)
    # Remove duplicate centroids (if necessary)
    # This can happen if the substructure matches overlap
    unique_centroids = []
    for centroid in match_centroids:
        if not any(np.allclose(centroid, uc) for uc in unique_centroids):
            unique_centroids.append(centroid)
    # Compute vectors from cage centroid to match centroids
    vectors = [centroid - cage_centroid for centroid in unique_centroids]
    distances = [np.linalg.norm(vector) for vector in vectors]
    unit_vectors = [
        vector / np.linalg.norm(vector)
        if np.linalg.norm(vector) != 0
        else vector
        for vector in vectors
    ]
    unit_vectors = np.array(unit_vectors)
    return unit_vectors, np.mean(distances)


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
    molecule: stk.Molecule,
    vectors: np.ndarray,
    vertice_size: float,
    no_vectors_define_facet: int,
    tolerance: float = 0.1, #The tolerance is a parameter that sets the acceptable margin of variation in the distances between vertices in a given combination – including slight irregularities – when determining whether they form a valid facet. For example if the distance between two vertices on a facet is 5A, and the tolerance is 0.1, then any vector which is less than 5.5A is considered to make up the facet. 
) -> np.ndarray:
    """Get axes using the midpoint of given vectors."""
    if isinstance(vectors, list) and isinstance(vectors[0], np.ndarray):
        vectors = vectors[
            0
        ]  # Assuming the first element is the NumPy array with all vectors
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

def remove_common(
    arr1: np.ndarray,
    arr2: np.ndarray,
    tolerance: float = 0.1,
) -> np.ndarray:
    """Remove common vectors."""
    filtered_arr1 = []
    for vec1 in arr1:
        diff = True
        for vec2 in arr2:
            similarity = cosine_similarity([vec1], [vec2])[0][0]
            if similarity > (1 - tolerance):
                diff = False
                break
        if diff:
            filtered_arr1.append(vec1)

    return np.array(filtered_arr1)


def remove_close_vectors(vectors, threshold=0.1):
    """Removes vectors from the list if the distance between any two vectors
    is less than the specified threshold.

    Parameters
    ----------
    vectors (numpy array): Array of vectors.
    threshold (float): Distance threshold below which one of the vectors will
    be removed.

    Returns
    -------
    numpy array: Array of vectors with close vectors removed.

    """
    # Convert vectors to a list to easily remove elements
    filtered_vectors = vectors.tolist()

    # Iterate over all vector pairs
    i = 0
    while i < len(filtered_vectors) - 1:
        j = i + 1
        while j < len(filtered_vectors):
            # Calculate the distance between the ith and jth vectors
            dist = np.linalg.norm(
                np.array(filtered_vectors[i]) - np.array(filtered_vectors[j])
            )
            if dist < threshold:
                # Remove the jth vector if it's too close to the ith vector
                filtered_vectors.pop(j)
            else:
                j += 1
        i += 1

    # Convert back to a numpy array before returning
    return np.array(filtered_vectors)
