import os
import stk
import stko
from . import utils
from . import config
from .cage import CageOperations
from .pores import *
import uuid
import pywindow as pw
import shutil
from itertools import combinations
from stko._internal.optimizers.utilities import (
    get_metal_atoms,
)
import rdkit.Chem.AllChem as rdkit
from sklearn.metrics.pairwise import cosine_similarity

GULP_PATH = config.GULP_PATH

class Axes:
    def ByPywindow(self,filename): #This doesnt work, not sure I understand python
        molsys = pw.MolecularSystem.load_file(filename)
        mol = molsys.system_to_molecule()
        windows=mol.calculate_windows()
        com=mol.calculate_centre_of_mass()
        adjusted_windows = windows - com
        return adjusted_windows
        #print("Performing ByPywindow calculation or operation.")
        #return "Result of ByPywindow"
    def BySmarts(self,smarts_string):
        print("Performing BySmarts calculation or operation.")
        return "Result of BySmarts"
        rdkit_mol = self.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_mol)
        centroid_smiles = [self.get_centroid(atom_ids=atom_ids) for atom_ids in rdkit_mol.GetSubstructMatches(query=rdkit.MolFromSmarts(smarts_string))]
        centroid_smiles = np.asarray(centroid_smiles)
        centroid_mol = self.get_centroid()
        distances = [np.linalg.norm(smile - centroid_mol) for smile in centroid_smiles]
        vectors = np.array([(smile - centroid_mol) / np.linalg.norm(smile - centroid_mol) for smile in centroid_smiles])
        return vectors,np.mean(distances)
    def BySmiles(self, smiles_string):

        print("Performing BySmiles calculation or operation.")

        rdkit_mol = self.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_mol)
        centroid_smiles = [self.get_centroid(atom_ids=atom_ids) for atom_ids in rdkit_mol.GetSubstructMatches(query=rdkit.MolFromSmiles(smiles_string))]
        centroid_smiles = np.asarray(centroid_smiles)
        centroid_mol = self.get_centroid()
        distances = [np.linalg.norm(smile - centroid_mol) for smile in centroid_smiles]
        vectors = np.array([(smile - centroid_mol) / np.linalg.norm(smile - centroid_mol) for smile in centroid_smiles])
        return vectors,np.mean(distances)

    def ByMidpoint(self,vectors,vertice_size,no_vectors_define_facet, tolerance=0.1):
        #if isinstance(vectors, np.ndarray):
        #    vectors = [vectors[i] for i in range(vectors.shape[0])]
        #print()
        if isinstance(vectors, list) and isinstance(vectors[0], np.ndarray):
            vectors = vectors[0]  # Assuming the first element is the NumPy array with all vectors

    # Convert numpy array of vectors into a list of numpy arrays
        vectors_list = [vectors[i] for i in range(vectors.shape[0])]
        all_distances = [utils.distance(v1, v2) for v1, v2 in combinations(vectors_list, 2)]
        min_distance = min(filter(lambda d: d > 0, all_distances))

        midpoints = []  # List to store all midpoints that meet the condition
        midpoint_size = []  # List to store all midpoints that meet the condition

        # Check each combination of n vectors
        for combination in combinations(vectors_list, no_vectors_define_facet):
            distances = [utils.distance(combination[i], combination[(i + 1) % no_vectors_define_facet]) for i in range(no_vectors_define_facet)]
            if max(distances) - min(distances) <= tolerance * min_distance:
                midpoint = np.mean(combination, axis=0)
                midpoint_size.append(np.linalg.norm(midpoint))
                midpoint=utils.normalize_vector(midpoint)
                midpoints.append(midpoint)  # Add the computed midpoint to the list

        return np.array(midpoints),np.mean(midpoint_size)*vertice_size

    def RemoveCommon(self,arr1, arr2,tolerance=0.1):

        filtered_arr1=[]
        for vec1 in arr1:
            diff = True
            for i, vec2 in enumerate(arr2):
                similarity = cosine_similarity([vec1], [vec2])[0][0]
                if similarity > (1-tolerance):
                    diff = False
                    break
            if diff:
                filtered_arr1.append(vec1)

        return np.array(filtered_arr1)


class DimerGenerator:

    def __init__(
        self,
        axes: np.ndarray,
        displacement: float = 7,
        displacement_step_size: float = 1,
        rotation_limit: float = 120,
        rotation_step_size: float = 30,
        overlap_tolerance: float = 0.1,
        slide: bool = False,
    ):
        self._axes = axes
        self._displacement = displacement
        self._displacement_step_size = displacement_step_size
        self._rotation_limit = rotation_limit
        self._rotation_step_size = rotation_step_size
        self._overlap_tolerance = overlap_tolerance
        self._slide = slide

    def generate(self,
        axes: np.ndarray,
        second_cage_orientation: np.ndarray,
        displacement_distance: float,
        displacement: float = 7,
        displacement_step_size: float = 1,
        rotation_limit: float = 120,
        rotation_step_size: float = 30,
        overlap_tolerance: float = 0.2,
        slide: bool = False):

        cage = stk.BuildingBlock.init_from_molecule(self)
        origin = cage.get_centroid()
        guest_cage=cage.with_rotation_between_vectors(second_cage_orientation,axes, origin)
        rotated_vectors=utils.generate_rotated_vectors(axes, rotation_limit/rotation_step_size, 30)
        perpendicular_vector=utils.calculate_perpendicular_vector(axes)

        dimer_list = []
        for i in range(0, int(displacement/displacement_step_size)):
            if slide:
                displaced_centers=utils.find_integer_points(axes, displacement_distance*2-2+i, int(displacement_distance) + 1)
            else:
                displaced_centers= [(displacement_distance*2-2+i)*axes]
            for center in displaced_centers:
                rot_by=0
                for vector in rotated_vectors:
                    rotated_guest = utils.create_rotated_guest(guest_cage,perpendicular_vector,vector,center)
                    dimer = stk.ConstructedMolecule(
                        topology_graph=stk.host_guest.Complex(host=self, guests=rotated_guest)
                    )
                    mol = dimer.to_rdkit_mol()
                    overlaps = utils.check_overlaps(mol,overlap_tolerance)
                    if overlaps:
                        rot_by=rot_by+1
                        continue
                    dimer_list.append({
                        'Displacement shell':(2-2+i),
                        'Displacement centroid': center,
                        'Rotation': rot_by*rotation_step_size,
                        'Dimer': dimer
                    })
                    rot_by=rot_by+1
        return dimer_list

        #def test_overlap


class Dimer(stko.GulpUFFOptimizer):

    def __init__(
        self,
        gulp_path: str,
        maxcyc: int = 1000,
        metal_FF: dict | None = None,
        metal_ligand_bond_order: str | None = None,
        conjugate_gradient: bool = False,
        output_dir: str | None = None,
        fixed_atom_set: list[int] | None = None,
    ):
        super().__init__(gulp_path=gulp_path, output_dir=output_dir, metal_FF=metal_FF,
                         conjugate_gradient=conjugate_gradient, maxcyc=maxcyc,metal_ligand_bond_order=metal_ligand_bond_order)
        self._fixed_atom_set = fixed_atom_set

    def optimize(self, mol: stk.Molecule,fixed_atom_set=None) -> stk.Molecule:
        if self._output_dir is None:
            output_dir = str(uuid.uuid4().int)
        else:
            output_dir = self._output_dir
        output_dir = os.path.abspath(output_dir)

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        os.mkdir(output_dir)
        init_dir = os.getcwd()
        os.chdir(output_dir)

        in_file = "gulp_opt.gin"
        out_file = "gulp_opt.ginout"
        output_xyz = "gulp_opt.xyz"

        metal_atoms = get_metal_atoms(mol)

        try:
            # Write GULP file.
            self._write_gulp_file(
                mol=mol,
                metal_atoms=metal_atoms,
                in_file=in_file,
                output_xyz=output_xyz,
                fixed_atom_set=fixed_atom_set,
            )
            # Run.
            self._run_gulp(in_file, out_file)

            # Update from output.
            mol = mol.with_structure_from_file(output_xyz)

        finally:
            os.chdir(init_dir)

        result=super().optimize(mol)

        return result

    def _constrain_section(self, fixed_atom_set) -> str:
        constrain_section = "\n"
        for constrain in fixed_atom_set:
            constrain_section += f"fix_atom {constrain}\n"

        return constrain_section
    def _write_gulp_file(
        self,
        mol: stk.Molecule,
        metal_atoms: list[stk.Atom],
        in_file: str,
        output_xyz: str,
        fixed_atom_set: list[int] | None = None,
    ) -> None:
        type_translator = self._type_translator()

        top_line = "opti "

        if self._conjugate_gradient:
            top_line += "conj unit "

        top_line += "conv "
        cell_section = ""
        periodic_output = ""

        top_line += "noautobond fix molmec cartesian\n"

        position_section = self._position_section(mol, type_translator)
        bond_section = self._bond_section(mol, metal_atoms)
        species_section = self._species_section(type_translator)
        constrain_section = self._constrain_section(fixed_atom_set)

        library = "\nlibrary uff4mof.lib\n"

        output_section = (
            "\n"
            f"maxcyc {self._maxcyc}\n"
            "terse inout potentials\n"
            "terse in cell\n"
            "terse in structure\n"
            "terse inout derivatives\n"
            f"output xyz {output_xyz}\n"
            f"{periodic_output}"
            # 'output movie xyz steps_.xyz\n'
        )

        with open(in_file, "w") as f:
            f.write(top_line)
            f.write(cell_section)
            f.write(position_section)
            f.write(bond_section)
            f.write(constrain_section)
            f.write(species_section)
            f.write(library)
            f.write(output_section)


class GulpDimerOptimizer:
    def __init__(self, gulp_path):
        self.gulp_path = gulp_path

    def optimise_dimer(self,cage, num, mode, dis_cent, rot, dis,fixed_atom_set):

        output_dir=f"Cage{num}_gulp/Cage{num}_{mode}/Cage_{num}_{dis_cent}_{rot}_{dis}_{mode}"
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, "gulp_opt.ginout")

        # Check if the file exists
        if os.path.exists(output_file):
            with open(output_file, 'r') as file:
                lines = file.readlines()
                # Check if the last non-empty line starts with the specified pattern
                for last_line in lines[::-1]:
                    if last_line.strip():  # This ensures we skip any empty lines at the end of the file
                        break

                if last_line.startswith("  Job Finished at "):
                    print(f"Skipping dimer Cage_{num}_{dis_cent}_{rot}_{dis}_{mode} as it is already done")
                    return  # Exit the function if the job is already done

        fix_atom = fixed_atom_set
        gulp_opt = Dimer(
            gulp_path=self.gulp_path,
            output_dir=output_dir,  # Change to correct path for Tmp files
            metal_FF={45: 'Rh6+3'},
            conjugate_gradient=True,
            maxcyc=500,
            fixed_atom_set=fix_atom,
        )
        gulp_opt.assign_FF(cage)
        structure = gulp_opt.optimize(mol=cage, fixed_atom_set=fix_atom)
        structure.write(f'{output_dir}_opt.mol')


# %%
    def run_calculations(self,Cagenum, molecule_type, molecule_list,fixed_atom_set):
        for molecule_cent in molecule_list:
            for molecule_rot in molecule_cent:
                for molecule in molecule_rot:
                    mol = molecule.to_rdkit_mol()
                    overlaps = utils.check_overlaps(mol)
                    if overlaps:
                        print(f"skip Cage_{Cagenum}_{molecule_cent.index(molecule_rot)}_{molecule_rot.index(molecule)}_{molecule_list.index(molecule_cent)}_{molecule_type}")
                        continue
                    utils.write_molecule_to_mol_file(molecule, Cagenum, molecule_type, molecule_cent.index(molecule_rot),molecule_rot.index(molecule),molecule_list.index(molecule_cent))
                    cage1=stk.BuildingBlock.init_from_file(f'Cage{Cagenum}/Cage{Cagenum}_{molecule_type}/Cage_{Cagenum}_{molecule_cent.index(molecule_rot)}_{molecule_rot.index(molecule)}_{molecule_list.index(molecule_cent)}_{molecule_type}.mol')
                    self.optimise_dimer(cage1,Cagenum, molecule_type,molecule_cent.index(molecule_rot),molecule_rot.index(molecule),molecule_list.index(molecule_cent),fixed_atom_set)


# %%
    def run_a_cage(self,Cagenum, cage, arene_smile, diamine_smile,sym="4_6",metal_atom=None,constrain=None,parallel=1,opt=False):
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

        ww=True
        wa=True
        aa=True
        wv=True

        conditions = {'ww': ww, 'wa': wa, 'aa': aa, 'wv': wv}
        arene_smile=utils.remove_aldehyde(arene_smile)
        diamine_smile=utils.remove_aldehyde(diamine_smile)

        for suffix in conditions:
            utils.create_folders_and_return_paths(foldername, [f"_{suffix}"])

        list_wa, list_ww, list_aa, list_wv, fixed_atom_set=CageOperations.displace_cages(cage, arene_smile, diamine_smile, sym,metal_atom)
        cage_lists = {'ww': list_ww, 'wa': list_wa, 'aa': list_aa, 'wv': list_wv}


        for suffix in cage_lists:
            specific_cage_list = cage_lists[suffix]

            # Convert the last cage of the specific list to an RDKit molecule
            rdkit_mol = specific_cage_list[-1][-1][-1].to_rdkit_mol()

            # Perform operations on the RDKit molecule
            fixed_atom_set = CageOperations.fix_atom_set(rdkit_mol, diamine_smile, metal_atom=metal_atom)

            # Run calculations for the specific cage list
            self.run_calculations(Cagenum, suffix, specific_cage_list, fixed_atom_set)

    def run_unconstrained_calculations(self,dir):
        for filename in os.listdir(dir):
        # Check if the file is a .mol file

            if filename.endswith(".mol"):
                output_dir=f"{dir}/{filename[:-4]}"
                os.makedirs(output_dir, exist_ok=True)
                output_file = os.path.join(output_dir, "gulp_opt.ginout")

        # Check if the file exists
                if os.path.exists(output_file):
                    with open(output_file, 'r') as file:
                        lines = file.readlines()
                        # Check if the last non-empty line starts with the specified pattern
                        for last_line in lines[::-1]:
                            if last_line.strip():  # This ensures we skip any empty lines at the end of the file
                                break

                        if last_line.startswith("  Job Finished at "):
                            print(f"Skipping dimer {filename[:-4]} as it is already done")
                        else: # Exit the function if the job is already done

                            # The path to your .mol file
                            mol_file_path = os.path.join(dir, filename)
                            cage=stk.BuildingBlock.init_from_file(mol_file_path)

                            gulp_opt = stko.GulpUFFOptimizer(
                                gulp_path=self.gulp_path,
                                output_dir=output_dir,  # Change to correct path for Tmp files
                                metal_FF={45: 'Rh6+3'},
                                conjugate_gradient=True,
                                maxcyc=500,
                            )
                            gulp_opt.assign_FF(cage)
                            structure = gulp_opt.optimize(mol=cage)
                            structure.write(f'{output_dir}_opt.mol')
                else:
                    mol_file_path = os.path.join(dir, filename)
                    cage=stk.BuildingBlock.init_from_file(mol_file_path)

                    gulp_opt = stko.GulpUFFOptimizer(
                        gulp_path=self.gulp_path,
                        output_dir=output_dir,  # Change to correct path for Tmp files
                        metal_FF={45: 'Rh6+3'},
                        conjugate_gradient=True,
                        maxcyc=500,
                    )
                    gulp_opt.assign_FF(cage)
                    structure = gulp_opt.optimize(mol=cage)
                    structure.write(f'{output_dir}_opt.mol')
    # Handle 'wv' condition separately

        #rdkit_mol = list_wa[-1][-1].to_rdkit_mol()
        #fixed_atom_set = CageOperations.fix_atom_set(rdkit_mol, diamine_smile,metal_atom=metal_atom)
#
        #self.run_calculations(Cagenum,'ww', list_ww,fixed_atom_set)
        #self.run_calculations(Cagenum,'wa', list_wa,fixed_atom_set)
        #self.run_calculations(Cagenum,'aa', list_aa,fixed_atom_set)
        #self.run_calculations(Cagenum,'wv', list_wv,fixed_atom_set)


class OPLSDimerOptimizer:

    def __init__(self, SCHRODINGER_PATH):
        self.SCHRODINGER_PATH = SCHRODINGER_PATH



    def make_files(self,Cagenum, molecule_type, molecule_list,fixed_atom_set,filename,cx1,constrain,parallel,opt):
        print(filename)
        if (cx1==False):
            with open(filename, 'a') as file:
                file.write(f'cd Cage{Cagenum}_xtb_{molecule_type} \n')


        for molecule_cent in molecule_list:
            for molecule_rot in molecule_cent:
                for molecule in molecule_rot:
                    mol = molecule.to_rdkit_mol()
                    overlaps = utils.check_overlaps(mol)
                    if overlaps:
                        print(f"skip Cage_{Cagenum}_{molecule_cent.index(molecule_rot)}_{molecule_rot.index(molecule)}_{molecule_list.index(molecule_cent)}_{molecule_type}")
                        continue
                    self.add_line_to_job_script(filename, f"Cage_{Cagenum}_{molecule_cent.index(molecule_rot)}_{molecule_rot.index(molecule)}_{molecule_list.index(molecule_cent)}_{molecule_type}",constrain,parallel,opt)
                    #self.write_molecule_to_mol_file_xtb(molecule, Cagenum, molecule_type, molecule_group.index(molecule), molecule_list.index(molecule_group))
                    self.write_molecule_to_mol_file_xtb(molecule, Cagenum, molecule_type, molecule_cent.index(molecule_rot),molecule_rot.index(molecule),molecule_list.index(molecule_cent))
                    #cage1=stk.BuildingBlock.init_from_file(f'Cage{Cagenum}/Cage{Cagenum}_ww/Cage_{Cagenum}_{molecule_group.index(molecule)}_{molecule_list.index(molecule_group)}_{molecule_type}.mol')
        if (cx1==False):
            with open(filename, 'a') as file:
                file.write('cd .. \n')

    def write_structcat_files(self,Cagenum,folder_path, foldername, molecule_type, molecule_list):
        with open(os.path.join(folder_path, foldername + ".txt"), "w") as file:
            file.write(f"export SCHRODINGER={self.SCHRODINGER_PATH}\n")
            file.write(" $SCHRODINGER/utilities/structcat")
            for molecule_cent in molecule_list:
                for molecule_rot in molecule_cent:
                    for molecule in molecule_rot:
                        mol = molecule.to_rdkit_mol()
                        overlaps = utils.check_overlaps(mol)
                        if overlaps:
                            print(f"skip Cage_{Cagenum}_{molecule_cent.index(molecule_rot)}_{molecule_rot.index(molecule)}_{molecule_list.index(molecule_cent)}_{molecule_type}")
                            continue
                        utils.write_molecule_to_mol_file(molecule, Cagenum, molecule_type, molecule_cent.index(molecule_rot),molecule_rot.index(molecule),molecule_list.index(molecule_cent))
                        file.write(f" -imae Cage_{molecule_cent.index(molecule_rot)}_{molecule_rot.index(molecule)}_{molecule_list.index(molecule_cent)}_{molecule_type}.mae")
            file.write(f" -omae {foldername}_{molecule_type}_merged.mae")
    # %%
    def run_a_cage(self,Cagenum, cage, arene_smile, diamine_smile,sym="4_6",metal_atom=None,constrain=None,parallel=1,opt=False):
        """
        The main function for running a calculation
        Args:
        - Cagenum: number of the cage
        - cage: optimized constructed cage
        - arene_smile: smiles string of

        """
        arene_smile=utils.remove_aldehyde(arene_smile)
        diamine_smile=utils.remove_aldehyde(diamine_smile)
        foldername = f"Cage{Cagenum}"
        os.makedirs(foldername, exist_ok=True)
        suffixes = ["_ww", "_wa", "_aa", "_wv"]
        utils.create_folders_and_return_paths(foldername, suffixes)
        foldername_mae = foldername + "_mae"
        os.makedirs(foldername_mae, exist_ok=True)
        mae_folder_paths = utils.create_folders_and_return_paths(foldername_mae, suffixes)


        list_wa, list_ww, list_aa, list_wv, fixed_atom_set=CageOperations.displace_cages(cage, arene_smile, diamine_smile, sym,metal_atom)

        new_content = utils.generate_com_content(fixed_atom_set)

        self.write_com_file(new_content, os.path.join(mae_folder_paths[0], foldername+".com"), foldername+"_ww_merged")
        self.write_com_file(new_content, os.path.join(mae_folder_paths[1], foldername+".com"), foldername+"_wa_merged")
        self.write_com_file(new_content, os.path.join(mae_folder_paths[2], foldername+".com"), foldername+"_aa_merged")
        self.write_com_file(new_content, os.path.join(mae_folder_paths[3], foldername+".com"), foldername+"_wv_merged")

        self.generate_sh_with_cage_number(foldername, f"{foldername_mae}/run_a_cage.sh")
        self.write_structcat_files(Cagenum,mae_folder_paths[0], foldername, 'ww', list_ww)
        self.write_structcat_files(Cagenum,mae_folder_paths[1], foldername, 'wa', list_wa)
        self.write_structcat_files(Cagenum,mae_folder_paths[2], foldername, 'aa', list_aa)
        self.write_structcat_files(Cagenum,mae_folder_paths[3], foldername, 'wv', list_wv)

    def generate_mode_files_and_script(self,foldername, mode, cages, cage_num, arene_smile, diamine_smile):
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
    export SCHRODINGER={self.SCHRODINGER_PATH}
    $SCHRODINGER/utilities/structcat {''.join(script_lines)} -omae {foldername}_{mode}_merged.mae
    """
        script_path = os.path.join(foldername, foldername + "_mae", mode_mae_folder, f"run_{foldername}_{mode}.sh")
        with open(script_path, 'w') as script_file:
            script_file.write(script_content)


    def create_directory_structure(self,base_folder):
        """
        Create the required directory structure for running calculations.

        Parameters:
        - base_folder: The base folder name where subfolders will be created.
        """
        modes = ['ww', 'wa', 'aa']
        for mode in modes:
            os.makedirs(os.path.join(base_folder, f'{base_folder}_{mode}'), exist_ok=True)
            os.makedirs(os.path.join(base_folder + "_mae", f'{base_folder}_mae_{mode}'), exist_ok=True)


    def write_com_files(self,folder_name, mode, content):
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

    SCHRODINGER_PATH = config.SCHRODINGER_PATH
    # %%
    def generate_com_content(self,fix_atoms):
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

    def generate_sh_with_cage_number(self,cage_number, output_file_name):
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

    export SCHRODINGER={self.SCHRODINGER_PATH}

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
    def write_com_file(self,new_content, filepath, name):
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

class XTBDimerOptimizer:
    def __init__(self, xtb_path):
        self.xtb_path = xtb_path


# %%
    def write_molecule_to_mol_file_xtb(self,molecule, num, mode, dis_cent, rot,dis):
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
            path=f'Cage{num}_xtb/Cage{num}_xtb_{mode}/Cage_{num}_{dis_cent}_{rot}_{dis}_{mode}.mol'
    )

    def add_line_to_job_script(self,filename, input_string,constrain,parallel,opt):
        # Constructing the line to be added
        line_to_add = "xtb "
        if constrain:
            line_to_add += "--input xtb.inp "
        line_to_add += f"--parallel {parallel} --namespace {input_string} {input_string}.mol --iterations 1000 "
        if opt:
            line_to_add += "--opt crude "
        line_to_add += f">{input_string}.out\n"

        # Opening the file in append mode to add the line
        with open(filename, 'a') as file:
            file.write(line_to_add)

    def make_files(self,Cagenum, molecule_type, molecule_list,fixed_atom_set,filename,cx1,constrain,parallel,opt):
        print(filename)
        if (cx1==False):
            with open(filename, 'a') as file:
                file.write(f'cd Cage{Cagenum}_xtb_{molecule_type} \n')


        for molecule_cent in molecule_list:
            for molecule_rot in molecule_cent:
                for molecule in molecule_rot:
                    mol = molecule.to_rdkit_mol()
                    overlaps = utils.check_overlaps(mol)
                    if overlaps:
                        print(f"skip Cage_{Cagenum}_{molecule_cent.index(molecule_rot)}_{molecule_rot.index(molecule)}_{molecule_list.index(molecule_cent)}_{molecule_type}")
                        continue
                    self.add_line_to_job_script(filename, f"Cage_{Cagenum}_{molecule_cent.index(molecule_rot)}_{molecule_rot.index(molecule)}_{molecule_list.index(molecule_cent)}_{molecule_type}",constrain,parallel,opt)
                    #self.write_molecule_to_mol_file_xtb(molecule, Cagenum, molecule_type, molecule_group.index(molecule), molecule_list.index(molecule_group))
                    self.write_molecule_to_mol_file_xtb(molecule, Cagenum, molecule_type, molecule_cent.index(molecule_rot),molecule_rot.index(molecule),molecule_list.index(molecule_cent))
                    #cage1=stk.BuildingBlock.init_from_file(f'Cage{Cagenum}/Cage{Cagenum}_ww/Cage_{Cagenum}_{molecule_group.index(molecule)}_{molecule_list.index(molecule_group)}_{molecule_type}.mol')
        if (cx1==False):
            with open(filename, 'a') as file:
                file.write('cd .. \n')

# %%

    def write_constraint_file(self,new_content, filepath):
        """
        Args:
        - new_content: the content to be written into the .com file
        - filepath: filepath
        - name: name of the file

        Returns:
        none
        """
        with open(filepath, 'w') as file:
            file.write(new_content)

    def generate_constraint_file(self,fix_atoms):
        combinations_list = combinations(fix_atoms, 2)
        formatted_combinations = ["    distance: {}, {}, auto".format(x, y) for x, y in combinations_list]

        # Constructing the new content with $constrain at the start and $end at the end
        new_content = "$constrain\n" + "\n".join(formatted_combinations) + "\n$end"
        return new_content


    def run_a_cage(self,Cagenum, cage, arene_smile, diamine_smile,sym="4_6",cx1=False,metal_atom=None,constrain=None,parallel=1,opt=False):
        """
        The main function for running a calculation
        Args:
        - Cagenum: number of the cage
        - cage: optimized constructed cage
        - arene_smile: smiles string of

        """

        arene_smile=utils.remove_aldehyde(arene_smile)
        diamine_smile=utils.remove_aldehyde(diamine_smile)
        #Create a folder named Cagenum
        foldername = f"Cage{Cagenum}"
        os.makedirs(foldername, exist_ok=True)
        suffixes = ["_ww", "_wa", "_aa", "_wv"]
        utils.create_folders_and_return_paths(foldername, suffixes)
        foldername_xtb = foldername + "_xtb"
        #os.makedirs(foldername_gulp, exist_ok=True)
        xtb_folder_paths = utils.create_folders_and_return_paths(foldername_xtb, suffixes)

        list_wa, list_ww, list_aa, list_wv, fixed_atom_set=CageOperations.displace_cages(cage, arene_smile, diamine_smile, sym,metal_atom)
        if constrain:
            new_content = self.generate_constraint_file(fixed_atom_set)
            self.write_constraint_file(new_content, os.path.join(xtb_folder_paths[0], "xtb.inp"))
            self.write_constraint_file(new_content, os.path.join(xtb_folder_paths[1], "xtb.inp"))
            self.write_constraint_file(new_content, os.path.join(xtb_folder_paths[2], "xtb.inp"))
            self.write_constraint_file(new_content, os.path.join(xtb_folder_paths[3], "xtb.inp"))
        if (cx1==True):
            self.write_job_script(foldername+"_ww",xtb_folder_paths[0],"constrain",parallel)
            full_path = f"{xtb_folder_paths[0]}/constrain.sh"
            self.make_files(Cagenum,'ww', list_ww,fixed_atom_set,full_path,cx1,constrain,parallel,opt)
            self.write_job_script(foldername+"_wa",xtb_folder_paths[1],"constrain",parallel)
            full_path = f"{xtb_folder_paths[1]}/constrain.sh"
            self.make_files(Cagenum,'wa', list_wa,fixed_atom_set,full_path,cx1,constrain,parallel,opt)
            self.write_job_script(foldername+"_aa",xtb_folder_paths[2],"constrain",parallel)
            full_path = f"{xtb_folder_paths[2]}/constrain.sh"
            self.make_files(Cagenum,'aa', list_aa,fixed_atom_set,full_path,cx1,constrain,parallel,opt)
            self.write_job_script(foldername+"_wv",xtb_folder_paths[3],"constrain",parallel)
            full_path = f"{xtb_folder_paths[3]}/constrain.sh"
            self.make_files(Cagenum,'wv', list_wv,fixed_atom_set,full_path,cx1,constrain,parallel,opt)
        else:
            filename=os.path.join(f'Cage{Cagenum}_xtb/', "constrain.sh")
            open(filename, 'w+')
            self.make_files(Cagenum,'ww', list_ww,fixed_atom_set,filename,cx1,constrain,parallel,opt)
            self.make_files(Cagenum,'wa', list_wa,fixed_atom_set,filename,cx1,constrain,parallel,opt)
            self.make_files(Cagenum,'aa', list_aa,fixed_atom_set,filename,cx1,constrain,parallel,opt)
            self.make_files(Cagenum,'wv', list_wv,fixed_atom_set,filename,cx1,constrain,parallel,opt)


    def write_job_script(self,job_name, filepath, filename,parallel):
        script_content = f"""#!/bin/bash --login
#PBS -N {job_name}
#PBS -l select=1:ncpus={parallel*2}:mem={parallel*2}gb:avx2=true
#PBS -l walltime=8:00:00

# Load modules for any applications

#module load mpi
#module load gcc
# Change to the directory the job was submitted from

cd $PBS_O_WORKDIR

# Run program, using 'mpiexec' to start the job
# mpiexec automatically picks up the # of cores
# assigned to the job. No other flags are required
#  - note: don't use 'mpirun'

#module load intel-suite
#module load cp2k/6.1-avx2
#module load gcc/10.2.0
#module load cuda/10.1
#module load mpi
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS={parallel}
export OMP_STACKSIZE=2G
export MKL_NUM_THREADS={parallel}
module load anaconda3/personal
source /rds/general/user/ewolpert/home/anaconda3/etc/profile.d/conda.sh
conda activate dimer_calculations
"""

        full_path = f"{filepath}/{filename}.sh"

        with open(full_path, 'w') as file:
            file.write(script_content)

