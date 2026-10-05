"""OpenMM energies of dimers and of the monomers they are built from.

All energies are in kJ/mol, from an OpenFF force field in vacuum:

``E_sp``
    The dimer exactly as built, without geometry optimisation.
``E_min``
    The dimer after minimisation.
``E_int_sp = E_sp - E(cage 1 as given) - E(cage 2 as given)``
    The interaction between the two cages without optimisation.
``E_int_min = E_min - E(cage 1 in min. dimer) - E(cage 2 in min. dimer)``
    The same interaction, but after minimisation.
``E_bind = E_min - E_min(cage 1 alone) - E_min(cage 2 alone)``
    The binding energy: the energy released when two relaxed cages form
    the relaxed dimer. It is the best single number for comparing motifs.
``E_def = E_bind - E_int_min``
    The deformation energy: how much more internal energy the two cages
    have in the minimised dimer than when relaxed alone.

Molecules sharing a topology (the same atoms in the same order, with the
same bonds) share one OpenMM context, so scoring many dimers of one cage
pair assigns force-field parameters only once.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
from openff.interchange import Interchange
from openff.toolkit import ForceField, Molecule, Topology
from openff.units import Quantity
from openmm import openmm, unit
from rdkit import Chem

if TYPE_CHECKING:
    from collections.abc import Sequence

    import stk

logger = logging.getLogger(__name__)

DEFAULT_FORCEFIELD = "openff_unconstrained-2.3.0.offxml"
_N_CAGES = 2


@dataclass(frozen=True)
class OpenMMSettings:
    """Settings for :class:`OpenMMCalculator`.

    Attributes
    ----------
    forcefield:
        An OpenFF force field. Sage 2.3.0, the default, assigns partial
        charges with NAGL, which takes milliseconds per molecule.

    platform:
        The OpenMM platform, for example ``"CPU"`` or ``"CUDA"``.

    max_iterations:
        Upper limit on L-BFGS iterations per minimisation.

    tolerance:
        Minimisation stops once the root-mean-square force falls below
        this value, in kJ/mol/nm.

    partial_charge_method:
        ``None`` uses the force field's own charges. Any other OpenFF
        method, such as ``"am1bcc"`` or ``"gasteiger"``, is computed once
        per molecule and reused.

    allow_undefined_stereo:
        Accept molecules with unassigned stereocentres.

    threads:
        Threads used by the CPU platform.

    """

    forcefield: str = DEFAULT_FORCEFIELD
    platform: str = "CPU"
    max_iterations: int = 2000
    tolerance: float = 10.0
    partial_charge_method: str | None = None
    allow_undefined_stereo: bool = True
    threads: int = 1


@dataclass
class _Context:
    """An OpenMM context and the atom order of its topology."""

    context: Any
    #: ``order[i]`` is the original index of topology atom ``i``.
    order: list[int]


def _to_rdkit(molecule: stk.Molecule) -> Chem.Mol:
    mol = molecule.to_rdkit_mol()
    Chem.SanitizeMol(mol)
    return mol


def _fragments(mol: Chem.Mol) -> list[tuple[int, ...]]:
    """Return the atom indices of each molecule, by first atom."""
    return sorted(
        Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False), key=min
    )


def _topology_key(mol: Chem.Mol) -> tuple:
    """Return elements plus the bond table; equal keys share a context."""
    elements = tuple(atom.GetAtomicNum() for atom in mol.GetAtoms())
    bonds = tuple(
        sorted(
            (
                bond.GetBeginAtomIdx(),
                bond.GetEndAtomIdx(),
                int(bond.GetBondTypeAsDouble() * 2),
            )
            for bond in mol.GetBonds()
        )
    )
    return elements, bonds


def _submol_in_order(mol: Chem.Mol, atom_indices: Sequence[int]) -> Chem.Mol:
    """Copy the given atoms, in the given order, into a new molecule."""
    rw = Chem.RWMol()
    mapping: dict[int, int] = {}
    for original in atom_indices:
        atom = mol.GetAtomWithIdx(int(original))
        new_atom = Chem.Atom(atom.GetAtomicNum())
        new_atom.SetFormalCharge(atom.GetFormalCharge())
        new_atom.SetChiralTag(atom.GetChiralTag())
        new_atom.SetIsAromatic(atom.GetIsAromatic())
        new_atom.SetNoImplicit(True)  # noqa: FBT003
        new_atom.SetNumExplicitHs(atom.GetNumExplicitHs())
        mapping[int(original)] = rw.AddAtom(new_atom)

    for bond in mol.GetBonds():
        begin, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if begin in mapping and end in mapping:
            rw.AddBond(mapping[begin], mapping[end], bond.GetBondType())

    sub = rw.GetMol()
    conformer = Chem.Conformer(len(atom_indices))
    source = mol.GetConformer()
    for original in atom_indices:
        conformer.SetAtomPosition(
            mapping[int(original)], source.GetAtomPosition(int(original))
        )
    sub.AddConformer(conformer, assignId=True)
    Chem.SanitizeMol(sub)
    return sub


class OpenMMCalculator:
    """Single points and minimisations with an OpenFF force field."""

    def __init__(self, settings: OpenMMSettings | None = None) -> None:
        """Initialise the calculator.

        Parameters
        ----------
        settings:
            Force field and minimiser settings. Defaults to
            :class:`OpenMMSettings`.

        """
        self.settings = settings or OpenMMSettings()
        self._contexts: dict[tuple, _Context] = {}
        self._charges: dict[str, np.ndarray] = {}

    def single_point(self, molecule: stk.Molecule) -> float:
        """Return the energy of ``molecule`` exactly as given.

        Parameters
        ----------
        molecule:
            The molecule to evaluate. No atom is moved.

        Returns
        -------
        float
            The potential energy, in kJ/mol.

        """
        return self._energy(_to_rdkit(molecule))

    def minimise(self, molecule: stk.Molecule) -> tuple[float, stk.Molecule]:
        """Minimise ``molecule`` with L-BFGS.

        Parameters
        ----------
        molecule:
            The starting geometry.

        Returns
        -------
        tuple[float, stk.Molecule]
            The minimised energy, in kJ/mol, and a copy of ``molecule``
            with the minimised geometry.

        """
        mol = _to_rdkit(molecule)
        held = self._load(mol)
        openmm.LocalEnergyMinimizer.minimize(
            held.context,
            tolerance=(
                self.settings.tolerance
                * unit.kilojoule_per_mole
                / unit.nanometer
            ),
            maxIterations=self.settings.max_iterations,
        )
        state = held.context.getState(getEnergy=True, getPositions=True)
        positions = np.empty((mol.GetNumAtoms(), 3))
        positions[held.order] = state.getPositions(asNumpy=True).value_in_unit(
            unit.angstrom
        )
        energy = state.getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole
        )
        return float(energy), molecule.with_position_matrix(positions)

    def fragment_energies(self, molecule: stk.Molecule) -> list[float]:
        """Return the energy of each molecule in ``molecule``, as placed.

        Parameters
        ----------
        molecule:
            For a dimer, the two cages.

        Returns
        -------
        list[float]
            One energy per molecule, in kJ/mol, ordered by first atom.

        """
        mol = _to_rdkit(molecule)
        return [
            self._energy(_submol_in_order(mol, fragment))
            for fragment in _fragments(mol)
        ]

    def _energy(self, mol: Chem.Mol) -> float:
        held = self._load(mol)
        state = held.context.getState(getEnergy=True)
        return float(
            state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        )

    def _load(self, mol: Chem.Mol) -> _Context:
        """Find or build the context for ``mol`` and set its positions."""
        key = _topology_key(mol)
        held = self._contexts.get(key)
        if held is None:
            held = self._build(mol)
            self._contexts[key] = held
        positions = mol.GetConformer().GetPositions()[held.order]
        held.context.setPositions(positions * unit.angstrom)
        return held

    def _build(self, mol: Chem.Mol) -> _Context:
        off_molecules = []
        preset: list[Molecule] = []
        order: list[int] = []
        for fragment in _fragments(mol):
            off_mol = Molecule.from_rdkit(
                _submol_in_order(mol, fragment),
                allow_undefined_stereo=self.settings.allow_undefined_stereo,
            )
            # Preset charges must be unique, so the second cage of a
            # homodimer is not listed again.
            if self._preset_charges(off_mol) and not any(
                off_mol.is_isomorphic_with(other) for other in preset
            ):
                preset.append(off_mol)
            off_molecules.append(off_mol)
            order.extend(int(index) for index in fragment)

        interchange = Interchange.from_smirnoff(
            force_field=ForceField(self.settings.forcefield),
            topology=Topology.from_molecules(off_molecules),
            charge_from_molecules=preset or None,
        )
        platform = openmm.Platform.getPlatformByName(self.settings.platform)
        properties = (
            {"Threads": str(self.settings.threads)}
            if self.settings.platform == "CPU"
            else {}
        )
        context = openmm.Context(
            interchange.to_openmm_system(),
            openmm.VerletIntegrator(1.0 * unit.femtoseconds),
            platform,
            properties,
        )
        logger.debug("Built an OpenMM context for %d atoms.", len(order))
        return _Context(context=context, order=order)

    def _preset_charges(self, off_mol: Molecule) -> bool:
        """Attach cached partial charges to ``off_mol``, if requested."""
        method = self.settings.partial_charge_method
        if method is None:
            return False
        key = off_mol.to_smiles(
            isomeric=True, explicit_hydrogens=True, mapped=True
        )
        if key not in self._charges:
            logger.info(
                "Assigning %s charges to a %d-atom molecule.",
                method,
                off_mol.n_atoms,
            )
            off_mol.assign_partial_charges(method)
            self._charges[key] = off_mol.partial_charges.m_as(
                "elementary_charge"
            )
        off_mol.partial_charges = Quantity(
            self._charges[key], "elementary_charge"
        )
        return True


@dataclass(frozen=True)
class MonomerEnergies:
    """A monomer's energy as given and after minimising it on its own.

    Attributes
    ----------
    as_given:
        Single-point energy of the monomer exactly as it is placed in the
        dimers, in kJ/mol.

    minimised:
        Energy after minimising the monomer on its own, in kJ/mol.

    """

    as_given: float
    minimised: float


def monomer_energies(
    monomer: stk.Molecule, calculator: OpenMMCalculator
) -> MonomerEnergies:
    """Calculate a monomer's reference energies.

    Parameters
    ----------
    monomer:
        The monomer exactly as it is used to build the dimers.

    calculator:
        The calculator, which is reused for the dimers.

    Returns
    -------
    MonomerEnergies
        The energy as given and after minimising on its own.

    """
    as_given = calculator.single_point(monomer)
    minimised, _ = calculator.minimise(monomer)
    return MonomerEnergies(as_given=as_given, minimised=minimised)


@dataclass(frozen=True)
class DimerEnergies:
    """The energies of one dimer, in kJ/mol.

    Attributes
    ----------
    e_sp:
        The dimer exactly as built, without geometry optimisation.

    e_min:
        The dimer after minimisation.

    cage_1_as_given, cage_2_as_given:
        Each monomer exactly as placed in the dimer.

    cage_1_in_min_dimer, cage_2_in_min_dimer:
        Each cage cut out of the minimised dimer, keeping that shape.

    cage_1_min_alone, cage_2_min_alone:
        Each monomer minimised on its own.

    minimised:
        The minimised dimer.

    """

    e_sp: float
    e_min: float
    cage_1_as_given: float
    cage_2_as_given: float
    cage_1_in_min_dimer: float
    cage_2_in_min_dimer: float
    cage_1_min_alone: float
    cage_2_min_alone: float
    minimised: stk.Molecule

    @property
    def e_int_sp(self) -> float:
        """Interaction between the cages without optimisation."""
        return self.e_sp - self.cage_1_as_given - self.cage_2_as_given

    @property
    def e_int_min(self) -> float:
        """Interaction between the cages after minimisation."""
        return self.e_min - self.cage_1_in_min_dimer - self.cage_2_in_min_dimer

    @property
    def e_bind(self) -> float:
        """Binding energy against monomers minimised on their own."""
        return self.e_min - self.cage_1_min_alone - self.cage_2_min_alone

    @property
    def e_def(self) -> float:
        """Deformation energy of both cages: ``e_bind - e_int_min``."""
        return self.e_bind - self.e_int_min

    def as_dict(self) -> dict[str, float]:
        """Return every energy under the usual table column name."""
        return {
            "E_bind": self.e_bind,
            "E_int_min": self.e_int_min,
            "E_int_sp": self.e_int_sp,
            "E_def": self.e_def,
            "E_min": self.e_min,
            "E_sp": self.e_sp,
            "E_cage1_in_min_dimer": self.cage_1_in_min_dimer,
            "E_cage2_in_min_dimer": self.cage_2_in_min_dimer,
            "E_cage1_as_given": self.cage_1_as_given,
            "E_cage2_as_given": self.cage_2_as_given,
            "E_min_cage1_alone": self.cage_1_min_alone,
            "E_min_cage2_alone": self.cage_2_min_alone,
        }


def dimer_energies(
    dimer: stk.Molecule,
    cage_1: MonomerEnergies,
    cage_2: MonomerEnergies,
    calculator: OpenMMCalculator,
) -> DimerEnergies:
    """Score one dimer as built, minimised, and against its monomers.

    The dimer must consist of the two monomers whose energies are given,
    each moved and turned rigidly, with the atoms of the first cage
    first. :class:`~dimer_calculations.dimer_generator.DimerGenerator`
    builds dimers this way.

    Parameters
    ----------
    dimer:
        The dimer as built.

    cage_1:
        Reference energies of the first monomer, from
        :func:`monomer_energies`.

    cage_2:
        Reference energies of the second monomer.

    calculator:
        The calculator used for the monomers.

    Returns
    -------
    DimerEnergies
        Every energy of the dimer, and its minimised geometry.

    """
    e_sp = calculator.single_point(dimer)
    e_min, minimised = calculator.minimise(dimer)
    in_dimer = calculator.fragment_energies(minimised)
    if len(in_dimer) != _N_CAGES:
        msg = f"Expected two cages in the dimer, found {len(in_dimer)}."
        raise ValueError(msg)
    return DimerEnergies(
        e_sp=e_sp,
        e_min=e_min,
        cage_1_as_given=cage_1.as_given,
        cage_2_as_given=cage_2.as_given,
        cage_1_in_min_dimer=in_dimer[0],
        cage_2_in_min_dimer=in_dimer[1],
        cage_1_min_alone=cage_1.minimised,
        cage_2_min_alone=cage_2.minimised,
        minimised=minimised,
    )
