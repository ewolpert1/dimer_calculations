"""Module for optimiser functions."""

import logging
import os
import shutil
import uuid
from itertools import combinations
from uuid import uuid4

import numpy as np
import stk
import stko
from stko._internal.optimizers.utilities import (
    get_metal_atoms,
    mol_from_mae_file,
    move_generated_macromodel_files,
)

from .utils import (
    calculate_perpendicular_vector,
    check_overlaps,
    create_rotated_guest,
    find_integer_points,
    generate_rotated_vectors,
)


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

    def generate(
        self,
        axes: np.ndarray,
        second_cage_orientation: np.ndarray,
        displacement_distance: float,
        displacement: float = 7,
        displacement_step_size: float = 1,
        rotation_limit: float = 120,
        rotation_step_size: float = 30,
        overlap_tolerance: float = 0.2,
        slide: bool = False,
    ):
        cage = stk.BuildingBlock.init_from_molecule(self)
        origin = cage.get_centroid()
        guest_cage = cage.with_rotation_between_vectors(
            second_cage_orientation, axes, origin
        )
        rotated_vectors = generate_rotated_vectors(
            axes, rotation_limit / rotation_step_size, 30
        )
        perpendicular_vector = calculate_perpendicular_vector(axes)

        dimer_list = []
        for i in range(int(displacement / displacement_step_size)):
            if slide:
                displaced_centers = find_integer_points(
                    axes,
                    displacement_distance * 2 - 2 + i,
                    int(displacement_distance) + 1,
                )
                print(slide)
            else:
                displaced_centers = [
                    (displacement_distance * 2 - 2 + i) * axes
                ]
            slide_up = 0
            for center in displaced_centers:
                rot_by = 0
                for vector in rotated_vectors:
                    rotated_guest = create_rotated_guest(
                        guest_cage, perpendicular_vector, vector, center
                    )
                    dimer = stk.ConstructedMolecule(
                        topology_graph=stk.host_guest.Complex(
                            host=self, guests=rotated_guest
                        )
                    )
                    mol = dimer.to_rdkit_mol()
                    overlaps = check_overlaps(mol, overlap_tolerance)

                    if overlaps:
                        rot_by = rot_by + 1
                        continue
                    dimer_list.append(
                        {
                            "Displacement shell": (2 - 2 + i),
                            "Slide": slide_up,
                            "Rotation": rot_by * rotation_step_size,
                            "Displacement centroid": center,
                            "Dimer": dimer,
                        }
                    )
                    # dimer_list.append(dimer_entry)
                    rot_by = rot_by + 1
                slide_up = slide_up + 1
        return dimer_list

        # def test_overlap


class OPLSDimer(stko.MacroModelForceField):
    def __init__(
        self,
        macromodel_path: str,
        output_dir: str | None = None,
        restricted: bool = False,
        timeout: float | None = None,
        force_field: int = 16,
        maximum_iterations: int = 2500,
        minimum_gradient: float = 0.05,
        fixed_atom_set: list[int] | None = None,
    ) -> None:
        # self._check_params(
        #    minimum_gradient=minimum_gradient,
        #    maximum_iterations=maximum_iterations,
        # )

        super().__init__(
            macromodel_path=macromodel_path,
            output_dir=output_dir,
            force_field=force_field,
            maximum_iterations=maximum_iterations,
            minimum_gradient=minimum_gradient,
            timeout=timeout,
        )
        self._fixed_atom_set = fixed_atom_set

    def _generate_com(
        self, mol: stk.Molecule, run_name: str, fixed_atom_set
    ) -> None:
        """Create a ``.com`` file for a MacroModel optimization.

        The created ``.com`` file fixes all bond parameters which were
        not added by :meth:`~.Topology.construct`. This means all bond
        distances, bond angles and torsional angles are fixed, except
        for cases where it involves a bond added by
        :meth:`.Topology.construct`.

        This fixing is implemented by creating a ``.com`` file with
        various "FX" commands written within its body.

        Parameters
        ----------
            mol:
                The molecule which is to be optimized.

            run_name:
                The name of the run. The files generated by this run will
                have this name.

        """
        # logger.debug(f'Creating .com file for "{mol}".')

        # This is the body of the ``.com`` file. The line that begins
        # and ends with exclamation lines is replaced with the various
        # commands that fix bond distances and angles.
        line1 = ("FFLD", self._force_field, 1, 0, 0, 1, 0, 0, 0)
        line2 = ("BGIN", 0, 0, 0, 0, 0, 0, 0, 0)
        line3 = ("READ", 0, 0, 0, 0, 0, 0, 0, 0)
        line4 = ("CONV", 2, 0, 0, 0, self._minimum_gradient, 0, 0, 0)
        line5 = ("MINI", 1, 0, self._maximum_iterations, 0, 0, 0, 0, 0)
        line6 = ("END", 0, 1, 0, 0, 0, 0, 0, 0)

        com_block = "\n".join(
            [
                self._get_com_line(*line1),
                self._get_com_line(*line2),
                self._get_com_line(*line3),
                "!!!BLOCK_OF_FIXED_PARAMETERS_COMES_HERE!!!",
                self._get_com_line(*line4),
                self._get_com_line(*line5),
                self._get_com_line(*line6),
            ]
        )

        # If `restricted` is ``False`` do not add a fix block.
        if fixed_atom_set == None:
            com_block = com_block.replace(
                "!!!BLOCK_OF_FIXED_PARAMETERS_COMES_HERE!!!\n", ""
            )
        else:
            # This function adds all the lines which fix bond distances
            # and angles into com_block.
            fix_block = ""
            for atom in fixed_atom_set:
                args = ("FXAT", atom, 0, 0, 0, 100, 0, 0, 0)
                fix_block += self._get_com_line(*args)
                fix_block += "\n"
            # fix_block = [
            # f" FXAT     {atom:>3}      0      0      0   100.0000     0.0000     0.0000     0.0000" for atom in fixed_atom_set
            # ]
            com_block = com_block.replace(
                "!!!BLOCK_OF_FIXED_PARAMETERS_COMES_HERE!!!\n", fix_block
            )

        # Writes the .com file.
        with open(f"{run_name}.com", "w") as com:
            # The first line holds the .mae file containing the
            # molecule to be optimized.
            com.write(f"{run_name}.mae\n")
            # The second line holds the name of the output file of the
            # optimization.
            com.write(f"{run_name}-out.maegz\n")
            # Next is the body of the .com file.
            com.write(com_block)

    def optimize(self, mol: stk.Molecule, fixed_atom_set) -> stk.Molecule:
        run_name = str(uuid4().int)
        if self._output_dir is None:
            output_dir = run_name
        else:
            output_dir = self._output_dir

        mol_path = f"{run_name}.mol"
        mae_path = f"{run_name}.mae"
        # First write a .mol file of the molecule.
        mol.write(mol_path)
        # MacroModel requires a ``.mae`` file as input.
        self._run_structconvert(mol_path, mae_path)
        # generate the ``.com`` file for the MacroModel run.
        self._generate_com(mol, run_name, fixed_atom_set)
        # Run the optimization.
        self._run_bmin(mol, run_name)
        # Get the ``.maegz`` optimization output to a ``.mae``.
        self._convert_maegz_to_mae(run_name)
        rdkit_opt_mol = mol_from_mae_file(mae_path)
        mol = mol.with_position_matrix(
            rdkit_opt_mol.GetConformer().GetPositions()
        )
        move_generated_macromodel_files(run_name, output_dir)
        return mol


class XTBDimer(stko.XTB):
    incomplete: set[stk.Molecule]

    def __init__(
        self,
        xtb_path: str,
        gfn_version: int = 2,
        output_dir: str | None = None,
        opt_level: str = "normal",
        max_runs: int = 2,
        calculate_hessian: bool = True,
        num_cores: int = 1,
        electronic_temperature: float = 300,
        solvent_model: str = "gbsa",
        solvent: str | None = None,
        solvent_grid: str = "normal",
        charge: int = 0,
        num_unpaired_electrons: int = 0,
        unlimited_memory: bool = False,
        fixed_atom_set: list[int] | None = None,
    ) -> None:
        super().__init__(
            xtb_path=xtb_path,
            gfn_version=gfn_version,
            output_dir=output_dir,
            opt_level=opt_level,
            max_runs=max_runs,
            calculate_hessian=calculate_hessian,
            num_cores=num_cores,
            electronic_temperature=electronic_temperature,
            solvent_model=solvent_model,
            solvent=solvent,
            solvent_grid=solvent_grid,
            charge=charge,
            num_unpaired_electrons=num_unpaired_electrons,
            unlimited_memory=unlimited_memory,
        )
        self._fixed_atom_set = fixed_atom_set

    def generate_constraint_file(self, fix_atoms):
        if fix_atoms:
            combinations_list = combinations(fix_atoms, 2)
            formatted_combinations = [
                f"    distance: {x}, {y}, auto" for x, y in combinations_list
            ]

            # Constructing the new content with $constrain at the start and $end at the end
            new_content = (
                "$constrain\n" + "\n".join(formatted_combinations) + "\n$end"
            )
        else:
            new_content = None
        return new_content

    def _write_detailed_control(self, fixed_atom_set) -> None:
        new_content = self.generate_constraint_file(fixed_atom_set)
        string = f"$gbsa\n   gbsagrid={self._solvent_grid}"

        with open("det_control.in", "w") as f:
            f.write(string)
            if new_content:
                f.write("\n")
                f.write(new_content)

    def _run_optimizations(
        self,
        mol: stk.Molecule,
        fixed_atom_set,
    ) -> tuple[stk.Molecule, bool]:
        """Run loop of optimizations on `mol` using xTB.

        Parameters
        ----------
            mol:
                The molecule to be optimized.

        Returns
        -------
            mol:
                The optimized molecule.

            opt_complete:
                Returns ``True`` if the calculation is complete and
                ``False`` if the calculation is incomplete.

        """
        for run in range(self._max_runs):
            xyz = f"input_structure_{run+1}.xyz"
            out_file = f"optimization_{run+1}.output"
            mol.write(xyz)
            self._write_detailed_control(fixed_atom_set)
            self._run_xtb(xyz=xyz, out_file=out_file)
            # Check if the optimization is complete.
            coord_file = "xtbhess.coord"
            output_xyz = "xtbopt.xyz"
            opt_complete = self._is_complete(out_file)
            if not opt_complete:
                if os.path.exists(coord_file):
                    # The calculation is incomplete.
                    # Update mol from xtbhess.coord and continue.
                    mol = mol.with_structure_from_file(coord_file)
                else:
                    # Update mol from xtbopt.xyz.
                    mol = mol.with_structure_from_file(output_xyz)
                    # If the negative frequencies are small, then GFN
                    # may not produce the restart file. If that is the
                    # case, exit optimization loop and warn.
                    self.incomplete.add(mol)
                    return mol, opt_complete
            else:
                # Optimization is complete.
                # Update mol from xtbopt.xyz.
                mol = mol.with_structure_from_file(output_xyz)
                break

        return mol, opt_complete

    def optimize(self, mol: stk.Molecule, fixed_atom_set) -> stk.Molecule:
        """Optimize `mol`.

        Parameters
        ----------
            mol:
                The molecule to be optimized.

        Returns
        -------
            mol:
                The optimized molecule.

        """
        # Remove mol from self.incomplete if present.
        if mol in self.incomplete:
            self.incomplete.remove(mol)

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

        try:
            mol, complete = self._run_optimizations(mol, fixed_atom_set)
        finally:
            os.chdir(init_dir)

        if not complete:
            self.incomplete.add(mol)
            logging.warning(f"Optimization is incomplete for {mol}.")

        return mol


class GulpDimer(stko.GulpUFFOptimizer):
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
        super().__init__(
            gulp_path=gulp_path,
            output_dir=output_dir,
            metal_FF=metal_FF,
            conjugate_gradient=conjugate_gradient,
            maxcyc=maxcyc,
            metal_ligand_bond_order=metal_ligand_bond_order,
        )
        self._fixed_atom_set = fixed_atom_set

    def optimize(self, mol: stk.Molecule, fixed_atom_set=None) -> stk.Molecule:
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

        return mol

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


class DimerOptimizer:
    def optimise_dimer_gulp(dimer, output_dir, gulp_path, fixed_atom_set=None):
        if fixed_atom_set is None:
            fixed_atom_set = []
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, "gulp_opt.ginout")
        if os.path.exists(output_file):
            with open(output_file) as file:
                lines = file.readlines()
                # Check if the last non-empty line starts with the specified pattern
                for last_line in lines[::-1]:
                    if last_line.strip():  # This ensures we skip any empty lines at the end of the file
                        break

                if last_line.startswith("  Job Finished at "):
                    print(f"Skipping dimer {output_dir} as it is already done")
                    return

        gulp_opt = GulpDimer(
            gulp_path=gulp_path,
            output_dir=output_dir,  # Change to correct path for Tmp files
            # metal_FF={45: 'Rh6+3'},
            conjugate_gradient=True,
            maxcyc=500,
            fixed_atom_set=fixed_atom_set,
        )
        gulp_opt.assign_FF(dimer)
        structure = gulp_opt.optimize(mol=dimer, fixed_atom_set=fixed_atom_set)
        structure.write(f"{output_dir}_opt.mol")

    def optimise_dimer_OPLS(
        dimer, output_dir, SCHRODINGER_PATH, fixed_atom_set=None
    ):
        os.makedirs(output_dir, exist_ok=True)
        # output_file = os.path.join(output_dir, "gulp_opt.ginout")
        # if os.path.exists(output_file):
        #    with open(output_file, 'r') as file:
        #        lines = file.readlines()
        #        # Check if the last non-empty line starts with the specified pattern
        #        for last_line in lines[::-1]:
        #            if last_line.strip():  # This ensures we skip any empty lines at the end of #the file
        #                break
        #
        #        if last_line.startswith("  Job Finished at "):
        #            print(f"Skipping dimer {output_dir} as it is already done")
        #            return

        OPLS_opt = OPLSDimer(
            macromodel_path=SCHRODINGER_PATH,
            output_dir=output_dir,  # Change to correct path for Tmp files
            # metal_FF={45: 'Rh6+3'},
            # conjugate_gradient=True,
            # maxcyc=500,
            fixed_atom_set=fixed_atom_set,
        )
        # OPLS_opt.assign_FF(dimer)
        structure = OPLS_opt.optimize(mol=dimer, fixed_atom_set=fixed_atom_set)
        structure.write(f"{output_dir}_opt.mol")

    def optimise_dimer_XTB(dimer, output_dir, XTB_PATH, fixed_atom_set=None):
        os.makedirs(output_dir, exist_ok=True)

        XTB_opt = XTBDimer(
            xtb_path=XTB_PATH,
            output_dir=output_dir,  # Change to correct path for Tmp files
            unlimited_memory=True,
            opt_level="crude",
            max_runs=1,
            calculate_hessian=True,
            fixed_atom_set=fixed_atom_set,
        )
        # OPLS_opt.assign_FF(dimer)
        structure = XTB_opt.optimize(mol=dimer, fixed_atom_set=fixed_atom_set)
        structure.write(f"{output_dir}_opt.mol")


# %%
