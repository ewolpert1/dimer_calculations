import numpy as np
import re
import shutil
import subprocess
import os
import stk
import sys
from .pores import *
from .dimer_centroids import *
from . import config
def find_last_segment_in_folder(folder_name):
    parts = folder_name.split('_')
    if len(parts) >= 3:
        return f"{parts[-1]}"
    ret
    
SCHRODINGER_PATH=config.SCHRODINGER_PATH

def read_energy(Cagename_mae, Cagename,current_directory,destination_folder_end):
    os.chdir(current_directory)
    destination_folder=f'{destination_folder_end}/{Cagename}'
    if not os.path.exists(destination_folder):
        os.makedirs(destination_folder)
    shutil.copy2("no_constraints.com", destination_folder)
    merge = open (f"{destination_folder}/merge.txt", "w+")
    #merge.write('export SCHRODINGER=/opt/schrodinger/suites2023-3/ \n$SCHRODINGER/utilities/structcat ')
    merge.write(f'export SCHRODINGER={SCHRODINGER_PATH} \n$SCHRODINGER/utilities/structcat ')
    all_minima = []  
    input1 = open("lowest_energies_50.txt", "w+")
    for dirpath, dirnames, filenames in os.walk("."):
        if dirpath == destination_folder:  # Skip destination folder
            continue
        folder_name = os.path.basename(dirpath)  # Extract folder name
        if folder_name.startswith(Cagename_mae):
            last_segment = find_last_segment_in_folder(folder_name)
            if (last_segment== None):
                continue
            for filename in filenames:
                if filename.endswith(".log"):  # Identify log files
                    minima_for_file = [] 
                    structure_energy_pairs = []
                    temp_energy = None
                    energy_counter = 0  # Counter to track every second 'Total Energy'
                    structure_name = None  
                    filepath = os.path.join(dirpath, filename)
                    output_filename = os.path.join(dirpath, "energies.txt")
                    with open(output_filename, "w+") as input1, open(filepath, 'r') as f:
                        lines = f.readlines() 
                        for i, line in enumerate(lines):
                            if "Structure name, if any, appears on next line:" in line:
                                i += 1  # Move to the next line to get the structure name
                                structure_name = lines[i].strip()
                                #print(structure_name)  # Print the structure name
                    
                            elif 'Total Energy' in line:
                                energy_counter += 1
                                if energy_counter % 2 == 0:  # Process only every second 'Total Energy'
                                    word = re.findall('\\d*\\.?\\d+', line)
                                    current_energy = float(word[0])
                                    if structure_name:
                                        #print(current_energy,structure_name)                                    
                                        minima_for_file.append((current_energy, structure_name))
                                        structure_name = None  # Reset structure_name for the next iteration
                                        temp_energy = None  # Reset temp_energy for the next iteration
                                        energy_counter = 0                         
                        for element in sorted(minima_for_file):
                            input1.write(str(element) + "\n")  
                    min_energy = min(minima_for_file, key=lambda x: x[0])[0]
                    sorted_minima_for_file = sorted(minima_for_file, key=lambda x: x[0])
                    reference_energy = None
                    for energy, structure_name in sorted_minima_for_file:
                        path, one_mol = generate_structure_path(Cagename, structure_name)
                        centroids = return_centroids(path)
                        cage_dimer = stk.BuildingBlock.init_from_file(path)
                        dimer = two_pores(cage_dimer, centroids[0], centroids[1])
                        cage = stk.BuildingBlock.init_from_file(f"{one_mol}.mol")
                        pore_cage = one_pore(cage)

                        if 'aa' in structure_name:
                            cutoff = 0
                        else:
                            cutoff = 0.2
                        catenated = check_catenane(pore_cage, dimer, cutoff)
                        print(catenated,one_mol) 

                        # Check if the structure is not catenated and set it as the reference if not already set
                        if not catenated:
                            if reference_energy is None:
                                reference_energy = energy  # Set the first uncatenated structure's energy as the reference
                            elif energy - reference_energy > 30:
                                break
                            #selected_structures.append((energy, structure_name))
                            if reference_energy is not None:
                        # If reference structure is found, check if the current structure is within 30 kJ/mol of the reference energy
                                print(catenated,one_mol)
                                mae_filename = f"{structure_name}.mae"
                                merge.write(f' -imae {mae_filename}')
                                src = os.path.join(dirpath, mae_filename)
                                dst = os.path.join(destination_folder, mae_filename)
                                shutil.copy2(src, dst)
                                all_minima.append((energy, structure_name))
                       
    merge.write(' -omae merged.mae \n')   
    os.chmod(f'{destination_folder}/merge.txt', os.stat(f'{destination_folder}/merge.txt').st_mode | 0o111)                
    return all_minima

def generate_structure_path(cage_number,structure_name):
    # Split the structure_name by underscores
    parts = structure_name.split('_')

    # Construct the first part using the first three segments
    first_part = '_'.join(parts[:4])

    # Construct the second part using the first six segments
    second_part =  (parts[-1])

    # Construct the final path
    path = f"{cage_number}/{cage_number}_{second_part}/{structure_name}.mol"
    one_mol=f'cages/{cage_number[4:]}'
    print(one_mol,path)    
    return path,one_mol

def run_calc(current_directory,Cagename,destination_folder_end):
    destination_folder=f'{destination_folder_end}/{Cagename}'    
    os.chdir(destination_folder)
    merge_command = "sh merge.txt"
    subprocess.run(merge_command, shell=True)    
    os.environ["SCHRODINGER"] = f"{SCHRODINGER_PATH}"
    schrodinger_bmin_command =  f"{SCHRODINGER_PATH}/bmin -WAIT no_constraints"
    subprocess.run(schrodinger_bmin_command, shell=True)  
    os.chdir(current_directory)

