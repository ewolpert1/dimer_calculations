"""Utilities module."""

import os
import pathlib
import subprocess

import numpy as np
import stk
from rdkit import Chem
from rdkit.Chem import rdMolTransforms


def load_matching_file(file_path):
    """Load and create a dictionary from a matching file."""
    mapping = {}
    with open(file_path) as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 2:
                mapping[parts[0]] = parts[1]
    return mapping


def remove_aldehyde(smiles_string: str) -> str:
    """Get the smiles string with aldehydes removed."""
    molecule = Chem.MolFromSmiles(smiles_string)
    aldehyde_smarts = "[C;H1](=O)"
    aldehyde_pattern = Chem.MolFromSmarts(aldehyde_smarts)
    updated_molecule = Chem.DeleteSubstructs(molecule, aldehyde_pattern)
    return Chem.MolToSmiles(updated_molecule, isomericSmiles=True)


def normalize_vector(vector: np.ndarray) -> np.ndarray:
    """Normalize a vector."""
    return vector / np.linalg.norm(vector)


def generate_com_content(fix_atoms):
    """Generate the content of a .com file used for minimization, including fixed atoms.

    Args:
    ----
    - fix_atoms: List of atom indices to be fixed.

    Returns:
    -------
    - A list of strings, each representing a line in the .com file.

    """
    header = [
        "merged.mae",
        " merged.maegz",
        " MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000",
        " DEBG      55      0      0      0     0.0000     0.0000     0.0000     0.0000",
        " FFLD      16      1      0      0     1.0000     0.0000     0.0000     0.0000",
        " BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000",
        " CRMS       0      0      0      0     0.0000     0.5000     0.0000     0.0000",
        " BGIN       0      0      0      0     0.0000     0.0000     0.0000     0.0000",
        " READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000",
    ]

    footer = [
        " CONV       2      0      0      0     0.0500     0.0000     0.0000     0.0000",
        " MINI       1      0   2500      0     0.0000     0.0000     0.0000     0.0000",
        " END        0      0      0      0     0.0000     0.0000     0.0000     0.0000",
    ]
    # Combine header, fxat lines for fixed atoms, and footer
    if fix_atoms:
        fxat_lines = [
            f" FXAT     {atom:>3}      0      0      0   100.0000     0.0000     0.0000     0.0000"
            for atom in fix_atoms
        ]
        com_content = header + fxat_lines + footer
    else:
        com_content = header + footer
    return com_content


def rotate_vector(v, axis, angle):
    """Rotate a vector around an axis by a given angle.

    Parameters
    ----------
    - v: Vector to be rotated.
    - axis: Axis of rotation.
    - angle: Angle of rotation in radians.

    Returns
    -------
    - Rotated vector.

    """
    axis = normalize_vector(axis)
    cos_angle = np.cos(angle)
    sin_angle = np.sin(angle)
    rotation_matrix = np.array(
        [
            [
                cos_angle + axis[0] ** 2 * (1 - cos_angle),
                axis[0] * axis[1] * (1 - cos_angle) - axis[2] * sin_angle,
                axis[0] * axis[2] * (1 - cos_angle) + axis[1] * sin_angle,
            ],
            [
                axis[1] * axis[0] * (1 - cos_angle) + axis[2] * sin_angle,
                cos_angle + axis[1] ** 2 * (1 - cos_angle),
                axis[1] * axis[2] * (1 - cos_angle) - axis[0] * sin_angle,
            ],
            [
                axis[2] * axis[0] * (1 - cos_angle) - axis[1] * sin_angle,
                axis[2] * axis[1] * (1 - cos_angle) + axis[0] * sin_angle,
                cos_angle + axis[2] ** 2 * (1 - cos_angle),
            ],
        ]
    )
    return np.dot(rotation_matrix, v)


def displace_vector(vector, distance):
    """Calculate a displacement vector along the given vector's direction with specified distance.

    Parameters
    ----------
    - vector: Direction vector for displacement.
    - distance: Magnitude of displacement.

    Returns
    -------
    - Displacement vector.

    """
    normalized_vector = normalize_vector(vector)
    return normalized_vector * distance


def find_integer_points(vector, dist, radius):
    """Find integer points within a circle defined along a vector.

    Parameters
    ----------
    - vector: Direction vector.
    - dist: Distance from the origin to the center of the circle.
    - radius: Radius of the circle.

    Returns
    -------
    - A list of integer points within the circle.

    """
    unit_vector = normalize_vector(vector)
    center = dist * unit_vector
    orthogonal_vectors = find_orthogonal_vectors(unit_vector)
    integer_points = []

    for x in range(-radius, radius + 1, int(radius / 2)):
        for y in range(-radius, radius + 1, int(radius / 2)):
            if np.linalg.norm([x, y]) <= radius:
                point3D = (
                    center
                    + x * orthogonal_vectors[0]
                    + y * orthogonal_vectors[1]
                )
                integer_points.append(point3D)

    return integer_points


def find_orthogonal_vectors(vector):
    """Find two orthogonal unit vectors perpendicular to a given vector.

    Parameters
    ----------
    - vector: The vector to find orthogonal vectors for.

    Returns
    -------
    - A tuple of two orthogonal unit vectors.

    """
    if vector[0] == 0 and vector[1] == 0:
        v1 = np.array([1, 0, 0])
    else:
        v1 = np.array([0, 0, 1])

    N1 = np.cross(vector, v1)
    N1 /= np.linalg.norm(N1)
    N2 = np.cross(vector, N1)
    N2 /= np.linalg.norm(N2)

    return N1, N2


def calculate_perpendicular_vector(vector):
    """Calculate a vector perpendicular to the given vector.

    Parameters
    ----------
    - vector: A 2D or 3D vector.

    Returns
    -------
    - A vector perpendicular to the input vector.

    """
    if len(vector) == 2:  # Two-dimensional case
        return np.array([-vector[1], vector[0]])
    elif len(vector) == 3:  # Three-dimensional case
        x, y, z = vector
        return np.array([y, -x, 0])
    else:
        raise ValueError("The input vector must be 2D or 3D.")


def create_rotated_guest(
    guest_building_block, start_vector, end_vector, displacement
):
    """Create a rotated and displaced guest molecule from a building block.

    Parameters
    ----------
    - guest_building_block: Building block from which the guest is created.
    - start_vector: Start vector for rotation.
    - end_vector: End vector for rotation.
    - displacement: Displacement vector from the building block.

    Returns
    -------
    - A guest molecule with specified rotation and displacement.

    """
    return stk.host_guest.Guest(
        building_block=guest_building_block,
        start_vector=start_vector,
        end_vector=end_vector,
        displacement=displacement,
    )


def distance(point1, point2):
    """Calculate the distance between two points."""
    return np.linalg.norm(point1 - point2)


def closest_point_on_segment(point, segment_start, segment_end):
    """Find the closest point on the line segment to the given point."""
    seg_vec = segment_end - segment_start
    pt_vec = point - segment_start
    seg_len = np.linalg.norm(seg_vec)
    seg_unit_vec = seg_vec / seg_len
    projection_length = np.dot(pt_vec, seg_unit_vec)
    if projection_length < 0:
        return segment_start
    if projection_length > seg_len:
        return segment_end
    return segment_start + projection_length * seg_unit_vec


def generate_rotated_vectors(base_vector, num_steps, angle_interval):
    """Generate a set of vectors rotated around the base vector.

    Parameters
    ----------
    - base_vector: The base vector around which the rotations will be performed.
    - num_steps: Number of vectors to generate.
    - angle_interval: Angle interval between each vector in degrees.

    Returns
    -------
    - An array of rotated vectors.

    """
    base_vector = normalize_vector(base_vector)
    perpendicular_vector = calculate_perpendicular_vector(base_vector)
    perpendicular_vector = normalize_vector(
        perpendicular_vector
        - np.dot(perpendicular_vector, base_vector) * base_vector
    )
    axis = base_vector
    angle_interval = np.radians(angle_interval)
    angles = np.arange(num_steps) * angle_interval
    return np.array(
        [rotate_vector(perpendicular_vector, axis, angle) for angle in angles]
    )


def create_folders_and_return_paths(parent_folder, suffixes):
    folder_paths = []  # Initialize a list to store the paths of the new folders
    for suffix in suffixes:
        # Construct the new folder path
        new_folder_path = os.path.join(parent_folder, parent_folder + suffix)
        # Create the new folder
        os.makedirs(new_folder_path, exist_ok=True)
        # Append the new folder path to the list
        folder_paths.append(new_folder_path)
    return folder_paths


def run_a_cage_script(cage_number):
    """Executes the run_a_cage.sh script for the specified cage number.

    Args:
    ----
    - cage_number: The specific cage number, used to construct the target directory path.

    """
    # Construct the target directory path
    current_dir = os.getcwd()
    target_directory = os.path.join(current_dir, f"Cage{cage_number}_mae")

    # Path to the run_a_cage.sh script
    script_path = os.path.join(
        target_directory, "run_a_cage.sh"
    )  # Adjust as necessary

    # Ensure the script is executable
    subprocess.run(["chmod", "+x", script_path], check=False)

    # Set the CAGE_DIRECTORY environment variable to the target directory
    env = os.environ.copy()
    env["CAGE_DIRECTORY"] = target_directory

    # Execute the shell script
    try:
        subprocess.run([script_path], check=True, env=env)
        print(f"Successfully executed run_a_cage.sh for {cage_number}")
    except subprocess.CalledProcessError as e:
        print(f"Error running run_a_cage.sh for {cage_number}: {e}")


def generate_rotated_vectors(base_vector, num_steps, angle_interval):
    """Generate a set of vectors rotated around the base vector.

    Parameters
    ----------
    - base_vector: The base vector around which the rotations will be performed.
    - num_steps: Number of vectors to generate.
    - angle_interval: Angle interval between each vector in degrees.

    Returns
    -------
    - An array of rotated vectors.

    """
    base_vector = normalize_vector(base_vector)
    perpendicular_vector = calculate_perpendicular_vector(base_vector)
    perpendicular_vector = normalize_vector(
        perpendicular_vector
        - np.dot(perpendicular_vector, base_vector) * base_vector
    )
    axis = base_vector
    angle_interval = np.radians(angle_interval)
    angles = np.arange(num_steps) * angle_interval
    return np.array(
        [rotate_vector(perpendicular_vector, axis, angle) for angle in angles]
    )


def write_molecule_to_mol_file(molecule, num, mode, dis_cent, rot, dis):
    """Save a given molecule to a .mol file with a specific naming convention.

    Parameters
    ----------
    - molecule: The stk molecule to be saved.
    - num: Identifier for the cage.
    - mode: The packing mode (e.g., 'wa', 'ww', 'aa').
    - dis: Displacement identifier.
    - rot: Rotation identifier.

    Example output file path: 'Cage100/Cage100_wa/Cage_100_1_1_wa.mol'

    """
    stk.MolWriter().write(
        molecule=molecule,
        path=f"Cage{num}/Cage{num}_{mode}/Cage_{num}_{dis_cent}_{rot}_{dis}_{mode}.mol",
    )


def mol_to_smiles(filepath: str | pathlib.Path) -> str:
    """Convert a .mol file to a SMILES string.

    Arguments:
    ---------
        filepath: Path to the .mol file.

    Returns:
    -------
    - A SMILES string representation of the molecule.

    """
    mol = Chem.MolFromMolFile(filepath)
    return Chem.MolToSmiles(mol)


def midpoint(conf: Chem.Mol, idx1: int, idx2: int) -> np.ndarray:
    """Calculate the midpoint of two atoms."""
    pos1 = np.array(conf.GetAtomPosition(idx1))
    pos2 = np.array(conf.GetAtomPosition(idx2))
    return (pos1 + pos2) / 2


def distance(point1, point2):
    """Calculate the distance between two points."""
    return np.linalg.norm(point1 - point2)


def closest_point_on_segment(point, segment_start, segment_end):
    """Find the closest point on the line segment to the given point."""
    seg_vec = segment_end - segment_start
    pt_vec = point - segment_start
    seg_len = np.linalg.norm(seg_vec)
    seg_unit_vec = seg_vec / seg_len
    projection_length = np.dot(pt_vec, seg_unit_vec)
    if projection_length < 0:
        return segment_start
    if projection_length > seg_len:
        return segment_end
    return segment_start + projection_length * seg_unit_vec


def check_overlaps(mol: Chem.Mol, threshold: float = 0.2) -> list[tuple]:
    """Check for overlaps between atoms in a molecule.

    Arguments:
    ---------
        mol (rdkit.Chem.Mol):
            The molecule to check.
        threshold (float):
            Distance threshold for overlap detection.

    Returns:
    -------
        overlaps (list):
            List of tuples with overlapping atom indices and their distance.

    """
    overlaps = []

    # Calculate pairwise distances
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        if mol.GetAtomWithIdx(i).GetAtomicNum() == 1:
            continue

        for j in range(i + 1, mol.GetNumAtoms()):
            # Skip hydrogen atoms
            if mol.GetAtomWithIdx(j).GetAtomicNum() == 1:
                continue

            dist = rdMolTransforms.GetBondLength(conf, i, j)

            if dist < threshold:
                overlaps.append((i, j, dist))
                # Remove this break if you want to find all overlaps for
                # each atom.
                break

    return overlaps
