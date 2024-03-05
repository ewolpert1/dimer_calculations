# utils.py
import os
from rdkit import Chem
import numpy as np
import stk
from rdkit import Chem
from rdkit.Chem import rdRGroupDecomposition, AllChem as rdkit, rdMolTransforms
import subprocess


def load_matching_file(file_path):
    """
    Load and create a dictionary from a matching file.
    """
    mapping = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 2:
                mapping[parts[0]] = parts[1]
    return mapping

def remove_aldehyde(smiles_string):
    """
    Given a trialdehyde smiles string, return the smiles string without aldehydes.
    """
    molecule = Chem.MolFromSmiles(smiles_string)
    aldehyde_smarts = '[C;H1](=O)'
    aldehyde_pattern = Chem.MolFromSmarts(aldehyde_smarts)
    updated_molecule = Chem.DeleteSubstructs(molecule, aldehyde_pattern)
    updated_smiles = Chem.MolToSmiles(updated_molecule, isomericSmiles=True)
    return updated_smiles

def normalize_vector(vector):
    """
    Normalize a vector.

    Parameters:
    - vector: A vector.

    Returns:
    - A normalized vector.
    """
    return vector / np.linalg.norm(vector)


def rotate_vector(v, axis, angle):
    """
    Rotate a vector around an axis by a given angle.

    Parameters:
    - v: Vector to be rotated.
    - axis: Axis of rotation.
    - angle: Angle of rotation in radians.

    Returns:
    - Rotated vector.
    """
    axis = normalize_vector(axis)
    cos_angle = np.cos(angle)
    sin_angle = np.sin(angle)
    rotation_matrix = np.array([
        [cos_angle + axis[0]**2 * (1 - cos_angle), 
         axis[0] * axis[1] * (1 - cos_angle) - axis[2] * sin_angle, 
         axis[0] * axis[2] * (1 - cos_angle) + axis[1] * sin_angle],
        [axis[1] * axis[0] * (1 - cos_angle) + axis[2] * sin_angle, 
         cos_angle + axis[1]**2 * (1 - cos_angle), 
         axis[1] * axis[2] * (1 - cos_angle) - axis[0] * sin_angle],
        [axis[2] * axis[0] * (1 - cos_angle) - axis[1] * sin_angle, 
         axis[2] * axis[1] * (1 - cos_angle) + axis[0] * sin_angle, 
         cos_angle + axis[2]**2 * (1 - cos_angle)]
    ])
    return np.dot(rotation_matrix, v)

def displace_vector(vector, distance):
    """
    Calculate a displacement vector along the given vector's direction with specified distance.

    Parameters:
    - vector: Direction vector for displacement.
    - distance: Magnitude of displacement.

    Returns:
    - Displacement vector.
    """
    normalized_vector = normalize_vector(vector)
    return normalized_vector * distance

def find_integer_points(vector, dist, radius):
    """
    Find integer points within a circle defined along a vector.

    Parameters:
    - vector: Direction vector.
    - dist: Distance from the origin to the center of the circle.
    - radius: Radius of the circle.

    Returns:
    - A list of integer points within the circle.
    """
    unit_vector = normalize_vector(vector)
    center = dist * unit_vector
    orthogonal_vectors = find_orthogonal_vectors(unit_vector)
    integer_points = []

    for x in range(-radius, radius + 1, int(radius / 2)):
        for y in range(-radius, radius + 1, int(radius / 2)):
            if np.linalg.norm([x, y]) <= radius:
                point3D = center + x * orthogonal_vectors[0] + y * orthogonal_vectors[1]
                integer_points.append(point3D)

    return integer_points

def find_orthogonal_vectors(vector):
    """
    Find two orthogonal unit vectors perpendicular to a given vector.

    Parameters:
    - vector: The vector to find orthogonal vectors for.

    Returns:
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
    """
    Calculate a vector perpendicular to the given vector.

    Parameters:
    - vector: A 2D or 3D vector.

    Returns:
    - A vector perpendicular to the input vector.
    """
    if len(vector) == 2:  # Two-dimensional case
        return np.array([-vector[1], vector[0]])
    elif len(vector) == 3:  # Three-dimensional case
        x, y, z = vector
        return np.array([y, -x, 0])
    else:
        raise ValueError("The input vector must be 2D or 3D.")
    
def create_rotated_guest(guest_building_block, start_vector, end_vector, displacement):
        """
        Create a rotated and displaced guest molecule from a building block.

        Parameters:
        - guest_building_block: Building block from which the guest is created.
        - start_vector: Start vector for rotation.
        - end_vector: End vector for rotation.
        - displacement: Displacement vector from the building block.

        Returns:
        - A guest molecule with specified rotation and displacement.
        """
        return stk.host_guest.Guest(
            building_block=guest_building_block,
            start_vector=start_vector,
            end_vector=end_vector,
            displacement=displacement,
        ) 
def midpoint(conf, idx1, idx2):
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

def check_overlaps(mol, threshold=0.8):
    """
    Check for overlaps between atoms in a molecule.
    
    Args:
    - mol (rdkit.Chem.Mol): The molecule to check.
    - threshold (float): Distance threshold for overlap detection.
    
    Returns:
    - overlaps (list): List of tuples with overlapping atom indices and their distance.
    """
    overlaps = []
    
    # Calculate pairwise distances
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        if mol.GetAtomWithIdx(i).GetAtomicNum() == 1:
            continue
            
        for j in range(i+1, mol.GetNumAtoms()):
            # Skip hydrogen atoms
            if mol.GetAtomWithIdx(j).GetAtomicNum() == 1:
                continue
            
            dist = rdMolTransforms.GetBondLength(conf, i, j)
            if dist < threshold:
                overlaps.append((i, j, dist))
                break  # Remove this break if you want to find all overlaps for each atom
            
                
    return overlaps

def generate_rotated_vectors(base_vector, num_steps, angle_interval):
    """
    Generate a set of vectors rotated around the base vector.

    Parameters:
    - base_vector: The base vector around which the rotations will be performed.
    - num_steps: Number of vectors to generate.
    - angle_interval: Angle interval between each vector in degrees.

    Returns:
    - An array of rotated vectors.
    """
    base_vector = normalize_vector(base_vector)
    perpendicular_vector = calculate_perpendicular_vector(base_vector)
    perpendicular_vector = normalize_vector(perpendicular_vector - np.dot(perpendicular_vector, base_vector) * base_vector)
    axis = base_vector
    angle_interval = np.radians(angle_interval)
    angles = np.arange(num_steps) * angle_interval
    return np.array([rotate_vector(perpendicular_vector, axis, angle) for angle in angles])

def run_a_cage_script(cage_number):
    """
    Executes the run_a_cage.sh script for the specified cage number.

    Args:
    - cage_number: The specific cage number, used to construct the target directory path.
    """

    # Construct the target directory path
    current_dir = os.getcwd()
    target_directory = os.path.join(current_dir, f"Cage{cage_number}_mae")

    # Path to the run_a_cage.sh script
    script_path = os.path.join(target_directory, "run_a_cage.sh")  # Adjust as necessary

    # Ensure the script is executable
    subprocess.run(["chmod", "+x", script_path])

    # Set the CAGE_DIRECTORY environment variable to the target directory
    env = os.environ.copy()
    env['CAGE_DIRECTORY'] = target_directory

    # Execute the shell script
    try:
        subprocess.run([script_path], check=True, env=env)
        print(f"Successfully executed run_a_cage.sh for {cage_number}")
    except subprocess.CalledProcessError as e:
        print(f"Error running run_a_cage.sh for {cage_number}: {e}")