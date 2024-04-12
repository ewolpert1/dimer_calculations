import os
import numpy as np
import pandas as pd
import stk
from rdkit import Chem
from rdkit.Chem import rdRGroupDecomposition, AllChem as rdkit, rdMolTransforms
import stko
from . import utils
from . import config
from .cage import Cage, CageOperations
from .pores import *
import re
import subprocess as sp
import uuid
import shutil
from itertools import combinations
from rdkit.Chem import AllChem as rdkit
from stko._internal.molecular.periodic.unitcell import UnitCell
from stko._internal.optimizers.optimizers import Optimizer
from stko._internal.optimizers.utilities import (
    get_metal_atoms,
    get_metal_bonds,
    has_h_atom,
    has_metal_atom,
    to_rdkit_mol_without_metals,
)
from stko._internal.utilities.exceptions import (
    ExpectedMetalError,
    ForceFieldSetupError,
    OptimizerError,
    PathError,
)

GULP_PATH = config.GULP_PATH


class GulpUFFOptimizer:
    """
    Applies forcefield optimizers that can handle metal centres.

    Notes:

        By default, :meth:`optimize` will run an optimisation using the
        UFF4MOF. This forcefield requires some explicit metal atom
        definitions, which are determined by the user.

        This code was originally written for use with Gulp 5.1 on Linux and has
        not been officially tested on other versions and operating systems.
        Make sure to sanity check the output.

    Examples:

        While metal atoms are not required, UFF4MOF is useful because it
        encompasses almost all chemical environments commonly found in
        metal-organic structures. Better forcefields exist for purely
        organic molecules! An interface with GULP is provided, which takes
        the forcefield types assigned by RDKit for non-metal atoms and
        user defined forcefield types for metal atoms to perform geometry
        optimisations.

        .. code-block:: python

            import stk
            import stko
            from rdkit.Chem import AllChem as rdkit

            # Produce a Pd+2 atom with 4 functional groups.
            atom = rdkit.MolFromSmiles('[Pd+2]')
            atom.AddConformer(rdkit.Conformer(atom.GetNumAtoms()))
            palladium_atom = stk.BuildingBlock.init_from_rdkit_mol(atom)
            atom_0, = palladium_atom.get_atoms(0)
            palladium_atom = palladium_atom.with_functional_groups(
                (stk.SingleAtom(atom_0) for i in range(4))
            )

            # Build a building block with two functional groups using
            # the SmartsFunctionalGroupFactory.
            bb1 = stk.BuildingBlock(
                smiles=('C1=CC(=CC(=C1)C2=CN=CC=C2)C3=CN=CC=C3'),
                functional_groups=[
                    stk.SmartsFunctionalGroupFactory(
                        smarts='[#6]~[#7X2]~[#6]',
                        bonders=(1, ),
                        deleters=(),
                    ),
                ],
            )

            # Build a metal-organic cage with dative bonds between
            # GenericFunctionalGroup and SingleAtom functional groups.
            cage = stk.ConstructedMolecule(
                stk.cage.M2L4Lantern(
                    building_blocks={
                        palladium_atom: (0, 1),
                        bb1: (2, 3, 4, 5)
                    },
                    reaction_factory=stk.DativeReactionFactory(
                        stk.GenericReactionFactory(
                            bond_orders={
                                frozenset({
                                    stk.GenericFunctionalGroup,
                                    stk.SingleAtom
                                }): 9
                            }
                        )
                    )
                )
            )

            # Perform Gulp optimisation with UFF4MOF.
            # Use conjugate gradient method for a slower, but more stable
            # optimisation.

            gulp_opt = stko.GulpUFFOptimizer(
                gulp_path='path/to/gulp',
                metal_FF={46: 'Pd4+2'},
                conjugate_gradient=True
            )

            # Assign the force field.
            gulp_opt.assign_FF(cage)
            # Run optimization.
            cage = gulp_opt.optimize(mol=cage)

    """

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
        """
        Parameters:

            gulp_path:
                Path to GULP executable.

            maxcyc:
                Set the maximum number of optimisation steps to use.
                Default in Gulp is 1000.

            metal_FF:
                Dictionary with metal atom forcefield assignments.
                Key: :class:`int` : atomic number.
                Value: :class:`str` : UFF4MOF forcefield type.

            metal_ligand_bond_order:
                Bond order to use for metal-ligand bonds. Defaults to
                `half`, but using `resonant` can increase the force
                constant for stronger metal-ligand interactions.

            conjugate_gradient:
                ``True`` to use Conjugate Graditent method.
                Defaults to ``False``

            output_dir:
                The name of the directory into which files generated during
                the calculation are written, if ``None`` then
                :func:`uuid.uuid4` is used.

        """

        self._check_path(gulp_path)
        self._gulp_path = gulp_path
        self._maxcyc = maxcyc
        self._metal_FF = metal_FF
        self._metal_ligand_bond_order = (
            "half"
            if metal_ligand_bond_order is None
            else metal_ligand_bond_order
        )
        self._conjugate_gradient = conjugate_gradient
        self._output_dir = output_dir
        self._fixed_atom_set = fixed_atom_set

    def _check_path(self, path: str) -> None:
        if not os.path.exists(path):
            raise PathError(f"GULP not found at {path}")

    def _add_atom_charge_flags(self, atom: rdkit.Atom, atomkey: str) -> str:
        """
        Add atom charge flags for forcefield.

        Code inspired by:
        https://github.com/rdkit/rdkit
        >   Code/GraphMol/ForceFieldHelpers/UFF/AtomTyper.cpp

        """
        total_valence = rdkit.Atom.GetTotalValence(atom)
        atnum = int(atom.GetAtomicNum())

        # Go through element cases.
        # Mg.
        if atnum == 12:
            if total_valence == 2:
                atomkey += "+2"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # Al.
        elif atnum == 13:
            if total_valence != 3:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )

        # Si.
        elif atnum == 14:
            if total_valence != 4:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # P.
        elif atnum == 15:
            if total_valence == 3:
                atomkey += "+3"
            elif total_valence == 5:
                atomkey += "+5"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )

        # S.
        elif atnum == 16:
            hybrid = rdkit.Atom.GetHybridization(atom)
            if hybrid != rdkit.HybridizationType.SP2:
                if total_valence == 2:
                    atomkey += "+2"
                elif total_valence == 4:
                    atomkey += "+4"
                elif total_valence == 6:
                    atomkey += "+6"
                else:
                    raise ForceFieldSetupError(
                        f"UFFTYPER: Unrecognized charge state for "
                        f"atom: {atom.GetIdx}"
                    )
        # Zn.
        elif atnum == 30:
            if total_valence == 2:
                atomkey += "+2"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )

        # Ga.
        elif atnum == 31:
            if total_valence == 3:
                atomkey += "+3"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # As.
        elif atnum == 33:
            if total_valence == 3:
                atomkey += "+3"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # Se.
        elif atnum == 34:
            if total_valence == 2:
                atomkey += "+2"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )

        # Cd.
        elif atnum == 48:
            if total_valence == 2:
                atomkey += "+2"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # In.
        elif atnum == 49:
            if total_valence == 3:
                atomkey += "+3"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )

        # Sb.
        elif atnum == 51:
            if total_valence == 3:
                atomkey += "+3"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # Te.
        elif atnum == 52:
            if total_valence == 2:
                atomkey += "+2"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # Hg.
        elif atnum == 80:
            if total_valence == 2:
                atomkey += "+2"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # Tl.
        elif atnum == 81:
            if total_valence == 3:
                atomkey += "+3"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # Pb.
        elif atnum == 82:
            if total_valence == 3:
                atomkey += "+3"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # Bi.
        elif atnum == 83:
            if total_valence == 3:
                atomkey += "+3"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # Po.
        elif atnum == 84:
            if total_valence == 2:
                atomkey += "+2"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # Lanthanides.
        elif atnum >= 57 and atnum <= 71:
            if total_valence == 6:
                atomkey += "+3"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        return atomkey

    def _get_atom_label(self, atom: rdkit.Atom) -> str:
        """
        Get FF atom label.

        Code inspired by:
        https://github.com/rdkit/rdkit
        >   Code/GraphMol/ForceFieldHelpers/UFF/AtomTyper.cpp

        """
        atnum = int(atom.GetAtomicNum())
        atomkey = atom.GetSymbol()
        if len(atomkey) == 1:
            atomkey += "_"

        table = rdkit.GetPeriodicTable()

        chk1 = rdkit.PeriodicTable.GetDefaultValence(table, atnum) == -1
        chk2 = rdkit.PeriodicTable.GetNOuterElecs(table, atnum) != 1
        chk3 = rdkit.PeriodicTable.GetNOuterElecs(table, atnum) != 7
        chk4 = chk2 and chk3
        if chk1 or chk4:
            hybrid = rdkit.Atom.GetHybridization(atom)
            if atnum == 84:
                atomkey += "3"
                if hybrid != rdkit.HybridizationType.SP3:
                    raise ForceFieldSetupError(
                        f"UFFTYPER: Unrecognized hybridization for"
                        f" atom: {atom.GetIdx}"
                    )
            elif atnum == 80:
                atomkey += "1"
                if hybrid != rdkit.HybridizationType.SP:
                    raise ForceFieldSetupError(
                        f"UFFTYPER: Unrecognized hybridization for"
                        f" atom: {atom.GetIdx}"
                    )
            else:
                if hybrid == rdkit.HybridizationType.SP:
                    atomkey += "1"
                elif hybrid == rdkit.HybridizationType.SP2:
                    chk1a = rdkit.Atom.GetIsAromatic(atom)
                    bonds = rdkit.Atom.GetBonds(atom)
                    conjugated = False
                    for bond in bonds:
                        if rdkit.Bond.GetIsConjugated(bond):
                            conjugated = True
                            break
                    chk2a = conjugated
                    chk3a = atnum in [6, 7, 8, 16]
                    chk4a = chk1a or chk2a
                    if chk4a and chk3a:
                        atomkey += "R"
                    else:
                        atomkey += "2"
                elif hybrid == rdkit.HybridizationType.SP3:
                    atomkey += "3"
                elif hybrid == rdkit.HybridizationType.SP3D:
                    atomkey += "5"
                elif hybrid == rdkit.HybridizationType.SP3D2:
                    atomkey += "6"
                else:
                    raise ForceFieldSetupError(
                        f"UFFTYPER: Unrecognized hybridization for"
                        f" atom: {atom.GetIdx}"
                    )
        atomkey = self._add_atom_charge_flags(atom, atomkey)
        return atomkey

    def _type_translator(self) -> dict[str, str]:
        type_translator: dict[str, str] = {}
        types = sorted(
            set(
                [  # type: ignore[type-var]
                    self.atom_labels[i][0] for i in self.atom_labels
                ]
            )
        )
        for t in types:
            if not t[1].isalpha():  # type:ignore[index]
                symb = t[0]  # type:ignore[index]
            else:
                symb = t[0:2]  # type:ignore[index]
            for i in range(1, 100):
                name = f"{symb}{i}"
                if name in type_translator.values():
                    continue
                else:
                    type_translator[t] = name  # type: ignore[index]
                    break

        return type_translator

    def _position_section(
        self, mol: stk.Molecule, type_translator: dict
    ) -> str:
        position_section = "\ncartesian\n"
        for atom in mol.get_atoms():
            atom_type = type_translator[self.atom_labels[atom.get_id()][0]]
            position = mol.get_centroid(atom_ids=atom.get_id())
            posi_string = (
                f"{atom_type} core {round(position[0], 5)} "
                f"{round(position[1], 5)} {round(position[2], 5)}\n"
            )
            position_section += posi_string

        return position_section


    def _bond_section(
        self,
        mol: stk.Molecule,
        metal_atoms: list[stk.Atom],
    ) -> str:
        bond_section = "\n"
        for bond in mol.get_bonds():
            atom_types = [
                self.atom_labels[i.get_id()][0]
                for i in [bond.get_atom1(), bond.get_atom2()]
            ]

            # Set bond orders.
            if has_h_atom(bond):
                # H has bond order of 1.
                bond_type = ""
            elif has_metal_atom(bond, metal_atoms):
                bond_type = self._metal_ligand_bond_order
            elif (
                "_R" in atom_types[0]  # type:ignore[operator]
                and "_R" in atom_types[1]  # type:ignore[operator]
            ):
                bond_type = "resonant"
            elif bond.get_order() == 1:
                bond_type = ""
            elif bond.get_order() == 2:
                bond_type = "double"
            elif bond.get_order() == 3:
                bond_type = "triple"

            string = (
                f"connect {bond.get_atom1().get_id()+1} "
                f"{bond.get_atom2().get_id()+1} {bond_type}"
            )
            bond_section += string + "\n"

        return bond_section

    def _species_section(self, type_translator: dict) -> str:
        species_section = "\nspecies\n"
        for spec in type_translator:
            name = type_translator[spec]
            species_section += f"{name} {spec}\n"

        return species_section

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

    def assign_FF(self, mol: stk.Molecule) -> None:
        """
        Assign forcefield types to molecule.

        Parameters:

            mol:
                The molecule to be optimized.

        """

        FutureWarning(
            "We have found some minor discrepancies in this "
            "assignment algorithm, which is based off rdkit code. "
            "Changes should come soon. This UFF optimisation should "
            " not be your final step! Due to this, some tests in "
            "test_uff_assign_ff.py have been muted."
        )

        metal_atoms = get_metal_atoms(mol)
        metal_ids = [i.get_id() for i in metal_atoms]

        if len(metal_ids) > 1 and self._metal_FF is None:
            raise ExpectedMetalError(
                "No metal FF provivded, but metal atoms were found ("
                f"{metal_atoms})"
            )

        metal_bonds, _ = get_metal_bonds(mol, metal_atoms)
        edit_mol = to_rdkit_mol_without_metals(
            mol=mol, metal_atoms=metal_atoms, metal_bonds=metal_bonds
        )

        # Get forcefield parameters.
        rdkit.SanitizeMol(edit_mol)
        self.atom_labels = {}

        for i in range(edit_mol.GetNumAtoms()):
            if i in metal_ids:
                self.atom_labels[i] = [None, "metal", None]
            else:
                atom = edit_mol.GetAtomWithIdx(i)
                atom_label = self._get_atom_label(atom)
                self.atom_labels[i] = [atom_label, None, None]

        # Write UFF4MOF specific forcefield parameters.
        # Metals.
        for atomid in self.atom_labels:
            if self.atom_labels[atomid][1] == "metal":
                (atom,) = mol.get_atoms(atomid)
                atom_no = atom.get_atomic_number()
                self.atom_labels[atomid][0] = (
                    self._metal_FF[  # type:ignore[index]
                        atom_no
                    ]
                )

    def _run_gulp(self, in_file: str, out_file: str) -> None:
        cmd = f"{self._gulp_path} < {in_file}"
        with open(out_file, "w") as f:
            # Note that sp.call will hold the program until completion
            # of the calculation.
            sp.call(
                cmd,
                stdin=sp.PIPE,
                stdout=f,
                stderr=sp.PIPE,
                # Shell is required to run complex arguments.
                shell=True,
            )


    def extract_final_energy(self, out_file: str) -> float:
        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        with open(out_file, "r") as f:
            for line in f.readlines():
                if "Final energy =" in line:
                    string = nums.search(line.rstrip())
                    return float(string.group(0))  # type: ignore[union-attr]

        raise OptimizerError(
            f'"Final energy =" not found in {out_file}, implying unsuccesful'
            " optimisation"
        )


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

        return mol


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
        gulp_opt = GulpUFFOptimizer(
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
        line_to_add = f"xtb "
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

