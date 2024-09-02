"""Read energies module."""

import glob
import os
import re
import shutil
import subprocess

import stk

from .config import SCHRODINGER_PATH
from .dimer_centroids import return_centroids
from .pores import check_catenane, one_pore, two_pores
from .utils import generate_com_content


def find_last_segment_in_folder(folder_name):
    parts = folder_name.split("_")
    if len(parts) >= 3:
        return f"{parts[-1]}"
    return None


class OPLSEnergyReader:
    def __init__(
        self,
        cage_name_mae,
        cage_name,
        current_directory,
        destination_folder_end,
    ):
        self.cage_name_mae = cage_name_mae
        self.cage_name = cage_name
        self.current_directory = current_directory
        self.destination_folder_end = destination_folder_end
        self.destination_folder = os.path.join(
            destination_folder_end, cage_name
        )

    def read_energy(self):
        os.chdir(self.current_directory)
        if not os.path.exists(self.destination_folder):
            os.makedirs(self.destination_folder)
        new_content = generate_com_content(fix_atoms=None)
        file_path = os.path.join(self.destination_folder, "no_constraints.com")
        with open(file_path, "w") as file:
            file.writelines("\n".join(new_content))
        merge_file_path = os.path.join(self.destination_folder, "merge.txt")
        with open(merge_file_path, "w+") as merge:
            merge.write(
                f"export SCHRODINGER={SCHRODINGER_PATH}\n$SCHRODINGER/utilities/structcat "
            )
        all_minima = []
        for dirpath, dirnames, filenames in os.walk(self.cage_name_mae):
            if dirpath == self.destination_folder:
                continue
            folder_name = os.path.basename(dirpath)  # Extract folder name
            if folder_name.startswith(self.cage_name_mae):
                last_segment = find_last_segment_in_folder(folder_name)
                if last_segment == None:
                    continue
                for filename in filenames:
                    if filename.endswith(".log"):  # Identify log files
                        self.process_log_file(
                            dirpath, filename, all_minima, merge_file_path
                        )
        self.finalize_merge_file(merge_file_path)
        sorted_minima = sorted(all_minima, key=lambda x: x[0])
        energies_file_path = os.path.join(
            self.destination_folder, "Constrained_energies_OPLS.txt"
        )
        with open(energies_file_path, "w") as energies_file:
            for energy, structure_name in sorted_minima:
                energies_file.write(f"{energy}, {structure_name}\n")
        return all_minima

    def process_log_file(self, dirpath, filename, all_minima, merge_file_path):
        minima_for_file = []
        energy_counter = 0
        structure_name = None
        filepath = os.path.join(dirpath, filename)
        output_filename = os.path.join(dirpath, "energies.txt")

        with open(output_filename, "w+") as input1, open(filepath) as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                if "Structure name, if any, appears on next line:" in line:
                    structure_name = lines[i + 1].strip()

                elif "Total Energy" in line:
                    energy_counter += 1
                    if energy_counter % 2 == 0:
                        energy_value = re.findall(r"\d*\.?\d+", line)[0]
                        current_energy = float(energy_value)
                        if structure_name:
                            minima_for_file.append(
                                (current_energy, structure_name)
                            )
                            structure_name = None
                            energy_counter = 0

            self.process_structures(
                minima_for_file, dirpath, merge_file_path, all_minima
            )

    def process_structures(self, minima, dirpath, merge_file_path, all_minima):
        sorted_minima = sorted(minima, key=lambda x: x[0])
        reference_energy = None

        for energy, structure_name in sorted_minima:
            path, one_mol = self.generate_structure_path(structure_name)
            centroids = return_centroids(path)
            cage_dimer = stk.BuildingBlock.init_from_file(path)
            dimer = two_pores(cage_dimer, centroids[0], centroids[1])
            cage = stk.BuildingBlock.init_from_file(f"{one_mol}.mol")
            pore_cage = one_pore(cage)

            cutoff = 0.1 if "aa" in structure_name else 0.2
            catenated = check_catenane(pore_cage, dimer, cutoff)

            if not catenated:
                if reference_energy is None:
                    # Set the first uncatenated structure's energy as the
                    # reference
                    reference_energy = energy
                elif energy - reference_energy > 30:
                    break
                # selected_structures.append((energy, structure_name))
                if reference_energy is not None:
                    self.update_merge_file(
                        merge_file_path, structure_name, dirpath
                    )
                    all_minima.append((energy, structure_name))

    def update_merge_file(self, merge_file_path, structure_name, dirpath):
        with open(merge_file_path, "a") as merge:
            mae_filename = f"{structure_name}.mae"
            merge.write(f" -imae {mae_filename}")
            src = os.path.join(dirpath, mae_filename)
            dst = os.path.join(self.destination_folder, mae_filename)
            shutil.copy2(src, dst)

    def finalize_merge_file(self, merge_file_path):
        with open(merge_file_path, "a") as merge:
            merge.write(" -omae merged.mae\n")
        os.chmod(merge_file_path, os.stat(merge_file_path).st_mode | 0o111)

    def generate_structure_path(self, structure_name):
        # Split the structure_name by underscores
        parts = structure_name.split("_")

        # Construct the second part using the first six segments
        second_part = parts[-1]

        # Construct the final path
        path = os.path.join(
            self.cage_name,
            f"{self.cage_name}_{second_part}",
            f"{structure_name}.mol",
        )
        one_mol = f"cages/{self.cage_name[4:]}"
        print(one_mol, path)
        return path, one_mol

    def run_calc(self):
        os.chdir(self.destination_folder)
        subprocess.run("sh merge.txt", shell=True, check=False)
        os.environ["SCHRODINGER"] = SCHRODINGER_PATH
        subprocess.run(
            f"{SCHRODINGER_PATH}/bmin -WAIT no_constraints",
            shell=True,
            check=False,
        )
        os.chdir(self.current_directory)


class GULPEnergyReader:
    def __init__(
        self,
        cage_name_gulp,
        cage_name,
        current_directory,
        destination_folder_end,
    ):
        self.cage_name_gulp = cage_name_gulp
        self.cage_name = cage_name
        self.current_directory = current_directory
        self.destination_folder_end = destination_folder_end
        self.destination_folder = os.path.join(
            destination_folder_end, cage_name
        )

    def read_energy(self, output_file="Constrained_energies_GULP.txt"):
        if not os.path.exists(self.destination_folder_end):
            os.makedirs(self.destination_folder_end)
            if not os.path.exists(self.destination_folder):
                os.makedirs(self.destination_folder)

        search_dir = self.cage_name_gulp
        all_minima = []
        for file_path in glob.glob(
            f"{search_dir}/**/gulp_opt.ginout", recursive=True
        ):
            folder_name = os.path.dirname(file_path)
            self.process_out_file(folder_name, file_path, all_minima)
        self.process_structures(all_minima)
        sorted_minima = sorted(all_minima, key=lambda x: x[0])
        energies_file_path = os.path.join(self.destination_folder, output_file)
        with open(energies_file_path, "w") as energies_file:
            for energy, structure_name in sorted_minima:
                energies_file.write(f"{energy}, {structure_name}\n")
        return all_minima
        #    if energy is not None:
        #        with open(output_file, 'a') as f:  # Open in append mode
        #            f.write(f"{folder_name}, {energy}\n")

    def process_out_file(self, folder_name, filename, all_minima):
        with open(filename) as f:
            for line in f:
                if "Final energy" in line:
                    parts = line.split()
                    if len(parts) >= 4:
                        if parts[3][0] != "*":
                            energy = float(
                                parts[3]
                            )  # Assuming the energy value is always in the 4th position
                            all_minima.append((energy, folder_name))
                    break

    def process_structures(self, minima):
        sorted_minima = sorted(minima, key=lambda x: x[0])
        reference_energy = None

        low_eng = []
        for energy, structure_name in sorted_minima:
            path, one_mol = self.generate_structure_path(structure_name)
            centroids = return_centroids(path)
            cage_dimer = stk.BuildingBlock.init_from_file(path)
            dimer = two_pores(cage_dimer, centroids[0], centroids[1])
            cage = stk.BuildingBlock.init_from_file(f"{one_mol}.mol")
            pore_cage = one_pore(cage)

            cutoff = 0.1 if "aa" in structure_name else 0.2
            catenated = check_catenane(pore_cage, dimer, cutoff)

            if not catenated:
                if reference_energy is None:
                    reference_energy = energy  # Set the first uncatenated structure's energy as the reference
                elif energy - reference_energy > 30:
                    break

                # selected_structures.append((energy, structure_name))
                if reference_energy is not None:
                    low_eng.append((energy, structure_name))
                    dst = os.path.join(
                        self.destination_folder, os.path.basename(path)
                    )
                    shutil.copy2(path, self.destination_folder)

    def generate_structure_path(self, structure_name):
        # Construct the final path
        path = f"{structure_name}_opt.mol"
        one_mol = f"cages/{self.cage_name[4:]}"
        print(one_mol, path)
        return path, one_mol


class XTBEnergyReader:
    def __init__(
        self,
        cage_name_xtb,
        cage_name,
        current_directory,
        destination_folder_end,
    ):
        self.cage_name_xtb = cage_name_xtb
        self.cage_name = cage_name
        self.current_directory = current_directory
        self.destination_folder_end = destination_folder_end
        self.destination_folder = os.path.join(
            destination_folder_end, cage_name
        )

    def read_energy(self, output_file="Constrained_energies_XTB.txt"):
        if not os.path.exists(self.destination_folder_end):
            os.makedirs(self.destination_folder_end)
            if not os.path.exists(self.destination_folder):
                os.makedirs(self.destination_folder)

        search_dir = self.cage_name_xtb
        all_minima = []
        for item in os.listdir(search_dir):
            print(item)
            dir_path = os.path.join(search_dir, item)

            if os.path.isdir(dir_path):
                print(dir_path)
                out_files = glob.glob(os.path.join(dir_path, "*.out"))
                for file_path in out_files:
                    folder_name = os.path.dirname(file_path)
                    self.process_out_file(folder_name, file_path, all_minima)
        self.process_structures(all_minima)
        sorted_minima = sorted(all_minima, key=lambda x: x[0])
        energies_file_path = os.path.join(self.destination_folder, output_file)
        with open(energies_file_path, "w") as energies_file:
            for energy, structure_name in sorted_minima:
                energies_file.write(f"{energy}, {structure_name}\n")
        return all_minima
        #    if energy is not None:
        #        with open(output_file, 'a') as f:  # Open in append mode
        #            f.write(f"{folder_name}, {energy}\n")

    def process_out_file(self, folder_name, filename, all_minima):
        with open(filename) as file:
            # Initialize a variable to hold the energy value
            energy = None

            # Loop through each line in the file
            for line in file:
                # Check if the line contains "TOTAL ENERGY"
                if "TOTAL ENERGY" in line:
                    energy = float(
                        line.split()[3]
                    )  # Assuming the energy value is the fourth element
                    # print(filename[:-4])
                    all_minima.append((energy, filename[:-4]))

            # If an energy value was found, write the filename and energy to the temp file

    def process_structures(self, minima):
        sorted_minima = sorted(minima, key=lambda x: x[0])
        reference_energy = None

        low_eng = []
        for energy, structure_name in sorted_minima:
            path, one_mol = self.generate_structure_path(structure_name)
            path_new = f"{structure_name}.mol"
            subprocess.run(["obabel", path, "-O", path_new, "-x3"], check=True)
            centroids = return_centroids(path_new)
            cage_dimer = stk.BuildingBlock.init_from_file(path_new)
            dimer = two_pores(cage_dimer, centroids[0], centroids[1])
            cage = stk.BuildingBlock.init_from_file(f"{one_mol}.mol")
            pore_cage = one_pore(cage)

            cutoff = 0.1 if "aa" in structure_name else 0.2
            catenated = check_catenane(pore_cage, dimer, cutoff)

            if not catenated:
                if reference_energy is None:
                    reference_energy = energy  # Set the first uncatenated structure's energy as the reference
                elif energy - reference_energy > 30:
                    break

                # selected_structures.append((energy, structure_name))
                if reference_energy is not None:
                    low_eng.append((energy, structure_name))
                    shutil.copy2(path_new, self.destination_folder)


    def generate_structure_path(self, structure_name):
        # Construct the final path
        path = f"{structure_name}.xtbopt.mol"
        one_mol = f"cages/{self.cage_name[4:]}"
        print(one_mol, path)
        return path, one_mol
