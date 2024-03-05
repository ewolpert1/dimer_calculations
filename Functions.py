import os
import numpy as np
import pandas as pd
import stk
from rdkit import Chem
from rdkit.Chem import rdRGroupDecomposition, AllChem as rdkit, rdMolTransforms
import stko
import utils
import config
from cage import Cage, CageOperations
from pores import *

# Set up Schrodinger path from configuration
SCHRODINGER_PATH = config.SCHRODINGER_PATH

# %%


# %%

# %%


# %%


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
    base_vector = utils.normalize_vector(base_vector)
    perpendicular_vector = utils.calculate_perpendicular_vector(base_vector)
    perpendicular_vector = utils.normalize_vector(perpendicular_vector - np.dot(perpendicular_vector, base_vector) * base_vector)
    axis = base_vector
    angle_interval = np.radians(angle_interval)
    angles = np.arange(num_steps) * angle_interval
    return np.array([utils.rotate_vector(perpendicular_vector, axis, angle) for angle in angles])




# %%

# %%

# %%
def write_molecule_to_mol_file(molecule, num, mode, dis, rot):
    """
    Save a given molecule to a .mol file with a specific naming convention.

    Parameters:
    - molecule: The stk molecule to be saved.
    - num: Identifier for the cage.
    - mode: The packing mode (e.g., 'wa', 'ww', 'aa').
    - dis: Displacement identifier.
    - rot: Rotation identifier.

    Example output file path: 'Cage100/Cage100_wa/Cage_100_1_1_wa.mol'
    """
    stk.MolWriter().write(
        molecule=molecule,
        path=f'Cage{num}/Cage{num}_{mode}/Cage_{num}_{dis}_{rot}_{mode}.mol'
    )


def cage_size(cage, arene_smile):
    """
    Calculate the average distance between the centroid of a cage and its arenes.

    Parameters:
    - cage: The stk cage molecule.
    - arene_smile: SMILES string representing the arene component.

    Returns:
    - The average distance between the centroid of the cage and its arenes.
    """
    rdkit_mol = cage.to_rdkit_mol()
    rdkit.SanitizeMol(rdkit_mol)
    arenes = [cage.get_centroid(atom_ids=atom_ids) for atom_ids in rdkit_mol.GetSubstructMatches(query=rdkit.MolFromSmarts(arene_smile))]
    centroid_mol = np.mean(arenes, axis=0)
    distances = [np.linalg.norm(arene - centroid_mol) for arene in arenes]
    return np.mean(distances)


# Creating a Cage object from a file


# Further operations can be called on the cage objects as needed



# %%
def fix_atom_set(rdkit_mol, diamine_smiles):
    """
    Return the atom IDs of the carbon atoms from the diamine that are adjacent to the amine nitrogen atoms.

    Args:
    - mol_file_path (str): Path to the .mol file.
    - diamine_smiles (str): SMILES string of the diamine. Default is "NC1CCCCC1N".

    Returns:
    - List[int]: List of atom IDs.
    """

    # Use the SMILES string of the diamine to identify its structure
    diamine = Chem.MolFromSmiles(diamine_smiles)

    # Get the substructure matches for the diamine
    matches = rdkit_mol.GetSubstructMatches(diamine)

    # List to store carbon atom ids
    carbon_ids = []

    # Iterate over the matches
    for match in matches:
        for atom_id in match:
            atom = rdkit_mol.GetAtomWithIdx(atom_id)
            # Check if the atom is carbon and adjacent to nitrogen
            if atom.GetSymbol() == "C" and any([neighbor.GetSymbol() == "N" for neighbor in atom.GetNeighbors()]):
                carbon_ids.append(atom_id)

    return carbon_ids


# %%
def generate_com_content(fix_atoms):
    """
    Generate the content of a .com file used for minimization, including fixed atoms.

    Args:
    - fix_atoms: List of atom indices to be fixed.

    Returns:
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

    fxat_lines = [
        f" FXAT     {atom:>3}      0      0      0   100.0000     0.0000     0.0000     0.0000" for atom in fix_atoms
    ]

    footer = [
        " CONV       2      0      0      0     0.0500     0.0000     0.0000     0.0000",
        " MINI       1      0   2500      0     0.0000     0.0000     0.0000     0.0000",
        " END        0      0      0      0     0.0000     0.0000     0.0000     0.0000"
    ]

    # Combine header, fxat lines for fixed atoms, and footer
    com_content = header + fxat_lines + footer

    return com_content


# %%

def generate_sh_with_cage_number(cage_number, output_file_name):
    """
    Generates a new .sh script for the specified cage number in the current directory.

    Args:
    - cage_number: Cage number e.g. "Cage100".
    - output_file_name: The name for the generated .sh file.

    Returns:
    None
    """
    current_dir = os.getcwd()
    script_content = f"""#!/bin/bash

# Source directory containing Mol files
source_dir_1="{current_dir}/{cage_number}/{cage_number}_ww"
source_dir_2="{current_dir}/{cage_number}/{cage_number}_wa"
source_dir_3="{current_dir}/{cage_number}/{cage_number}_aa"
source_dir_4="{current_dir}/{cage_number}/{cage_number}_wv"

# Destination directory to save MAE files
destination_dir_1="{current_dir}/{cage_number}_mae/{cage_number}_mae_ww"
destination_dir_2="{current_dir}/{cage_number}_mae/{cage_number}_mae_wa"
destination_dir_3="{current_dir}/{cage_number}_mae/{cage_number}_mae_aa"
destination_dir_4="{current_dir}/{cage_number}_mae/{cage_number}_mae_wv"

export SCHRODINGER={SCHRODINGER_PATH}

process_directory() {{
    source_dir=$1
    destination_dir=$2
    echo "Processing directory $source_dir"

    for mol_file in "$source_dir"/*.mol; do
        mol_filename=$(basename "$mol_file")
        mol_basename="${{mol_filename%.mol}}"
        mae_file="$destination_dir/$mol_basename.mae"

        if [ ! -f "$mae_file" ]; then
            $SCHRODINGER/utilities/structconvert "$mol_file" "$mae_file"
            echo "Converted $mol_filename to $mae_file"
        else
            echo "$mae_file already exists. Skipping conversion."
        fi
    done

    if [ ! "$(ls "$destination_dir"/*.log 2>/dev/null)" ]; then
        perform_tasks=true
    else
        perform_tasks=false
        for log_file in "$destination_dir"/*.log; do
            if ! grep -q "BatchMin: normal termination" "$log_file"; then
                perform_tasks=true
                break
            fi
        done
    fi

    if [ "$perform_tasks" = true ]; then
        cd "$destination_dir"
        python3 {current_dir}/Structure_name.py "$destination_dir"
        chmod +x {cage_number}.txt
        ./{cage_number}.txt
        $SCHRODINGER/bmin -WAIT {cage_number}
    fi
}}


# Process each directory
process_directory "$source_dir_1" "$destination_dir_1"
process_directory "$source_dir_2" "$destination_dir_2"
process_directory "$source_dir_3" "$destination_dir_3"
process_directory "$source_dir_4" "$destination_dir_4"
"""
    
    output_file_path = os.path.join(current_dir, output_file_name)
    # Write the script content to the output file
    with open(output_file_path, "w") as file:
        file.write(script_content)

    print(f"Generated script at {output_file_path} for {cage_number}.")


# %%
def write_com_file(new_content, filepath, name):
    """
    Args:
    - new_content: the content to be written into the .com file
    - filepath: filepath
    - name: name of the file

    Returns:
    none
    
    """
    with open(filepath, 'w') as file:
        file.write(name+".mae\n")
        file.write(name+".maegz\n")
        for line in new_content[2:]:
            file.write(line + "\n")

# %%
def mol_to_smiles(filepath):
    """
    Converts a .mol file to a SMILES string.
    
    Args:
    - filepath: Path to the .mol file.
    
    Returns:
    - A SMILES string representation of the molecule.
    """
    mol = Chem.MolFromMolFile(filepath)
    return Chem.MolToSmiles(mol)

# %%
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

# %%
def optimize_constructed_cage(bb1_smile, bb2_smile, aldehyde, amine):
    """
    Optimize the structure of a constructed cage molecule using a series of optimization steps.

    Parameters:
    - bb1_smile: SMILES string of the first building block.
    - bb2_smile: SMILES string of the second building block.
    - aldehyde: Identifier for the aldehyde component.
    - amine: Identifier for the amine component.

    Returns:
    - cage2: The optimized stk constructed cage molecule.
    - arene_smile: The SMILES string of the arene building block.
    - diamine_smile: The SMILES string of the diamine building block.
    """
    directories = ["cages", "MD", "opt", "FF_Restricted", "FF_Unrestricted"]
    for directory in directories:
        os.makedirs(directory, exist_ok=True)       
    cage_name=f"{aldehyde}_{amine}"
    cage_path = f'cages/{cage_name}.mol'

    arene_smile, diamine_smile = (bb1_smile, bb2_smile) if "C=O" in bb1_smile else (bb2_smile, bb1_smile)
    arene_smile = utils.remove_aldehyde(arene_smile)

    if os.path.exists(cage_path):
        cage2 = stk.BuildingBlock.init_from_file(cage_path)
    else:
        cage2 = construct_and_optimize_cage(bb1_smile, bb2_smile, cage_name)
        stk.MolWriter().write(cage2, cage_path)

    return cage2, arene_smile, diamine_smile


def construct_and_optimize_cage(bb1_smile, bb2_smile, cage_name):
    """
    Construct and optimize a cage molecule using a sequence of optimizers.

    Parameters:
    - bb1_smile: SMILES string of the first building block.
    - bb2_smile: SMILES string of the second building block.
    - cage_name: Name of the cage molecule for directory and file naming.

    Returns:
    - The optimized cage molecule.
    """
    cage = Cage()
    cage.construct(bb1_smile, bb2_smile)

    optimizer_sequence = stko.OptimizerSequence(
        stko.MacroModelForceField(f"FF_Restricted/{cage_name}", SCHRODINGER_PATH, 16, restricted=True),
        stko.MacroModelForceField(f"FF_Unrestricted/{cage_name}", SCHRODINGER_PATH, 16, restricted=False),
        stko.MacroModelMD(SCHRODINGER_PATH, f"MD/{cage_name}", temperature=700, conformers=50, simulation_time=100000, time_step=1, eq_time=100),
    )

    optimized_cage = optimizer_sequence.optimize(cage)
    return optimized_cage.with_centroid([0, 0, 0])
# %%
def run_a_cage(Cagenum, cage, arene_smile, diamine_smile):
    """
    The main function for running a calculation
    Args: 
    - Cagenum: number of the cage
    - cage: optimized constructed cage
    - arene_smile: smiles string of 
    
    """
    #Create a folder named Cagenum
    foldername = f"Cage{Cagenum}"
    os.makedirs(foldername, exist_ok=True)
    #Create folders inside Cagenum
    new_folder_path = os.path.join(foldername, foldername + "_ww")
    os.makedirs(new_folder_path, exist_ok=True)
    new_folder_path = os.path.join(foldername, foldername + "_wa")
    os.makedirs(new_folder_path, exist_ok=True)
    new_folder_path = os.path.join(foldername, foldername + "_aa")
    os.makedirs(new_folder_path, exist_ok=True)
    new_folder_path = os.path.join(foldername, foldername + "_wv")
    os.makedirs(new_folder_path, exist_ok=True)
    #Create a folder named Cagenum_mae
    foldername_mae = foldername+"_mae"
    os.makedirs(foldername_mae, exist_ok=True)
    #Create folders inside Cagenum_mae
    new_folder_path_1 = os.path.join(foldername_mae, foldername_mae+"_ww")
    os.makedirs(new_folder_path_1, exist_ok=True)
    new_folder_path_2 = os.path.join(foldername_mae, foldername_mae+"_wa")
    os.makedirs(new_folder_path_2, exist_ok=True)
    new_folder_path_3 = os.path.join(foldername_mae, foldername_mae+"_aa")
    os.makedirs(new_folder_path_3, exist_ok=True)
    new_folder_path_4 = os.path.join(foldername_mae, foldername_mae+"_wv")
    os.makedirs(new_folder_path_4, exist_ok=True)

    #generate all the runs for this specific cage

    vec = CageOperations.calculate_cage_vector(cage, arene_smile)
    dist = cage_size(cage, arene_smile)
    vec_vertex = cage_vertex_vec(cage, diamine_smile) #this has all the vectors of the arene to centroid
    guest_bb_wa = stk.BuildingBlock.init_from_molecule(cage)
    origin = guest_bb_wa.get_centroid()
    guest_bb_ww = guest_bb_wa.with_rotation_between_vectors(vec, vec*-1, origin)
    guest_bb_wv = guest_bb_wa.with_rotation_between_vectors(vec_vertex, vec, origin)
    #aa buildingblock
    guest_bb_aa = guest_bb_ww
    #vec_perpendicular
    vec_perpendicular = utils.calculate_perpendicular_vector(vec)
    #generate rotated vectors (rotate 4 times, angle invertal is 10 degrees)
    rotated_vectors = generate_rotated_vectors(vec, 4, 30)
    #generate all the runs for this specific cage
    list_wa, list_ww, list_aa, list_centroids = CageOperations.generate_new_cages(cage, guest_bb_aa, guest_bb_ww, guest_bb_wa, vec, vec_perpendicular, rotated_vectors, dist)
    list_wv, list_centroids = CageOperations.generate_ww_cages(guest_bb_wa, guest_bb_ww, vec, vec_perpendicular, rotated_vectors, dist)


    #generate a sh file in Cagenum_mae
    #write txt files and mol files


    #generate com file into each _mae folder
    rdkit_mol = list_wa[-1][-1].to_rdkit_mol()
    fixed_atom_set = fix_atom_set(rdkit_mol, diamine_smile)
    new_content = generate_com_content(fixed_atom_set)
    write_com_file(new_content, os.path.join(new_folder_path_1, foldername+".com"), foldername+"_ww_merged")
    write_com_file(new_content, os.path.join(new_folder_path_2, foldername+".com"), foldername+"_wa_merged")
    write_com_file(new_content, os.path.join(new_folder_path_3, foldername+".com"), foldername+"_aa_merged")
    write_com_file(new_content, os.path.join(new_folder_path_4, foldername+".com"), foldername+"_wv_merged")

    #generate a sh file in Cagenum_mae
    generate_sh_with_cage_number(foldername, f"{foldername_mae}/run_a_cage.sh")
    #write txt files and mol files
    with open(os.path.join(new_folder_path_1, foldername+".txt"), "w") as file:
        file.write(f"export SCHRODINGER={SCHRODINGER_PATH}"+"\n")
        file.write(" $SCHRODINGER/utilities/structcat")
        for ww in list_ww:
            for w_w in ww:
                mol = w_w.to_rdkit_mol()
                overlaps = check_overlaps(mol)
                if overlaps:
                    print(f"skip Cage_{Cagenum}_{ww.index(w_w)}_{list_ww.index(ww)}_ww")
                    continue
                write_molecule_to_mol_file(w_w, Cagenum, 'ww', ww.index(w_w), list_ww.index(ww))
                file.write(f" -imae Cage_{Cagenum}_{ww.index(w_w)}_{list_ww.index(ww)}_ww.mae")
        file.write(f" -omae {foldername}_ww_merged.mae")
    with open(os.path.join(new_folder_path_2, foldername+".txt"), "w") as file:
        file.write(f"export SCHRODINGER={SCHRODINGER_PATH}"+"\n")
        file.write(" $SCHRODINGER/utilities/structcat")
        for wa in list_wa:
            for w_a in wa:
                mol = w_a.to_rdkit_mol()
                overlaps = check_overlaps(mol)
                if overlaps:
                    print(f"skip Cage_{Cagenum}_{wa.index(w_a)}_{list_wa.index(wa)}_wa")
                    continue
                write_molecule_to_mol_file(w_a, Cagenum, 'wa', wa.index(w_a), list_wa.index(wa))
                file.write(f" -imae Cage_{Cagenum}_{wa.index(w_a)}_{list_wa.index(wa)}_wa.mae")
        file.write(f" -omae {foldername}_wa_merged.mae")
    with open(os.path.join(new_folder_path_3, foldername+".txt"), "w") as file:
        file.write(f"export SCHRODINGER={SCHRODINGER_PATH}"+"\n")
        file.write(" $SCHRODINGER/utilities/structcat")
        for aa in list_aa:
            for a_a in aa:
                mol = a_a.to_rdkit_mol()
                overlaps = check_overlaps(mol)
                if overlaps:
                    print(f"skip Cage_{Cagenum}_{aa.index(a_a)}_{list_aa.index(aa)}_aa")
                    continue
                write_molecule_to_mol_file(a_a, Cagenum, 'aa', aa.index(a_a), list_aa.index(aa))
                file.write(f" -imae Cage_{Cagenum}_{aa.index(a_a)}_{list_aa.index(aa)}_aa.mae")
        file.write(f" -omae {foldername}_aa_merged.mae")
    with open(os.path.join(new_folder_path_4, foldername+".txt"), "w") as file:
        file.write(f"export SCHRODINGER={SCHRODINGER_PATH}"+"\n")
        file.write(" $SCHRODINGER/utilities/structcat")
        for ww in list_wv:
            for w_w in ww:
                mol = w_w.to_rdkit_mol()
                overlaps = check_overlaps(mol)
                if overlaps:
                    print(f"skip Cage_{Cagenum}_{ww.index(w_w)}_{list_wv.index(ww)}_wv")
                    continue
                write_molecule_to_mol_file(w_w, Cagenum, 'wv', ww.index(w_w), list_wv.index(ww))
                file.write(f" -imae Cage_{Cagenum}_{ww.index(w_w)}_{list_wv.index(ww)}_wv.mae")
        file.write(f" -omae {foldername}_wv_merged.mae")        

def generate_mode_files_and_script(foldername, mode, cages, cage_num, arene_smile, diamine_smile):
    script_lines = []  # Collect lines for the script file
    mode_folder = f"{foldername}_{mode}"
    mode_mae_folder = f"{foldername}_mae_{mode}"

    for index, cage in enumerate(cages):
        mol_filename = f"Cage{cage_num}_{index}_{mode}.mol"
        mol_path = os.path.join(foldername, mode_folder, mol_filename)
        
        # Write the MOL file for the cage
        with open(mol_path, 'w') as mol_file:
            mol_file.write(cage.to_mol())

        # Prepare a line for the script to include this MOL file
        script_lines.append(f" -imae {mol_filename}")

    # Generate the script file for merging all MOL files into a single MAE
    script_content = f"""#!/bin/bash
export SCHRODINGER={SCHRODINGER_PATH}
$SCHRODINGER/utilities/structcat {''.join(script_lines)} -omae {foldername}_{mode}_merged.mae
"""
    script_path = os.path.join(foldername, foldername + "_mae", mode_mae_folder, f"run_{foldername}_{mode}.sh")
    with open(script_path, 'w') as script_file:
        script_file.write(script_content)


def create_directory_structure(base_folder):
    """
    Create the required directory structure for running calculations.

    Parameters:
    - base_folder: The base folder name where subfolders will be created.
    """
    modes = ['ww', 'wa', 'aa']
    for mode in modes:
        os.makedirs(os.path.join(base_folder, f'{base_folder}_{mode}'), exist_ok=True)
        os.makedirs(os.path.join(base_folder + "_mae", f'{base_folder}_mae_{mode}'), exist_ok=True)


def write_com_files(folder_name, mode, content):
    """
    Write .com files based on the provided content.

    Parameters:
    - folder_name: The base folder name where the .com file will be written.
    - mode: The packing mode ('ww', 'wa', 'aa').
    - content: The content of the .com file as a list of lines.
    """
    file_path = os.path.join(folder_name + "_mae", f'{folder_name}_mae_{mode}', f'{folder_name}_{mode}_merged.com')
    with open(file_path, 'w') as file:
        file.writelines('\n'.join(content))


# Note: The `MOL_writer` function might need adjustment to accept a file path parameter for writing the .mol file.



def cage_vertex_vec(stk_mol, amine_smile):
    """
    Returns the vector of average distance between the centroid of the cage and arenes

    """
    
    #pdb_file = "%s.pdb"%(cagename)
    rdkit_mol = stk_mol.to_rdkit_mol()
    rdkit.SanitizeMol(rdkit_mol)
    arenes = []
    for atom_ids in rdkit_mol.GetSubstructMatches(
        query=rdkit.MolFromSmarts(amine_smile),
    ):
        centroid = stk_mol.get_centroid(atom_ids=atom_ids)
        arenes.append(centroid)
    arenes=np.asarray(arenes)
    centroid_mol=stk_mol.get_centroid()
    vec=np.zeros((len(arenes),3))
    for i in range(len(arenes)):
        vec[i] = (arenes[i]-centroid_mol)
        vec[i]=vec[i]/(np.sqrt(np.dot(vec[i],vec[i])))
    return vec[0]





# %%
