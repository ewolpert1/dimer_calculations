"""Tests for OpenMM energies. Skipped when OpenMM or OpenFF is missing."""

import pathlib

import pytest
import stk

pytest.importorskip("openmm")
pytest.importorskip("openff.toolkit")

from dimer_calculations.axes import (
    by_midpoint,
    by_smiles,
    remove_common,
)
from dimer_calculations.dimer_generator import DimerGenerator
from dimer_calculations.openmm_energies import (
    OpenMMCalculator,
    OpenMMSettings,
    dimer_energies,
    monomer_energies,
)
from dimer_calculations.utils import remove_aldehyde

CC1 = pathlib.Path(__file__).parent.parent / "cages" / "CC1.mol"


@pytest.fixture(scope="module")
def cage() -> stk.BuildingBlock:
    """Return the CC1 cage, centred at the origin."""
    return stk.BuildingBlock.init_from_file(CC1).with_centroid([0, 0, 0])


@pytest.fixture(scope="module")
def calculator() -> OpenMMCalculator:
    """Return a calculator with a short minimisation, for quick tests."""
    return OpenMMCalculator(OpenMMSettings(max_iterations=200))


@pytest.fixture(scope="module")
def dimer(cage: stk.BuildingBlock) -> stk.Molecule:
    """Build one window-to-window CC1 dimer, as in the smoke test."""
    list_of_vertices, vertice_size = by_smiles(
        molecule=cage, smiles_string=remove_aldehyde("NCCN")
    )
    facet_axes, _centroid_to_facet = by_midpoint(
        molecule=cage,
        vectors=list_of_vertices,
        vertice_size=vertice_size,
        no_vectors_define_facet=3,
        tolerance=0.1,
    )
    list_of_arenes, centroid_to_arene = by_smiles(
        molecule=cage,
        smiles_string=remove_aldehyde("O=Cc1cc(C=O)cc(C=O)c1"),
    )
    list_of_windows = remove_common(facet_axes, list_of_arenes)

    generator = DimerGenerator(
        displacement=1,
        displacement_step_size=1,
        rotation_limit=30,
        rotation_step_size=30,
        slide=False,
    )
    dimers = generator.generate(
        molecule=cage,
        molecule_2=cage,
        axes_1=list_of_windows[0],
        axes_2=list_of_windows[0],
        displacement_distance=2 * centroid_to_arene + 2,
    )
    return dimers[0]["Dimer"]


def test_minimising_lowers_the_energy(
    cage: stk.BuildingBlock, calculator: OpenMMCalculator
) -> None:
    """Minimisation never raises the energy and keeps every atom."""
    start = calculator.single_point(cage)
    energy, minimised = calculator.minimise(cage)
    assert energy <= start
    assert minimised.get_num_atoms() == cage.get_num_atoms()
    # The returned geometry evaluates back to the minimised energy.
    assert calculator.single_point(minimised) == pytest.approx(
        energy, abs=1e-6
    )


def test_cages_of_a_built_dimer_have_the_monomer_energy(
    cage: stk.BuildingBlock,
    dimer: stk.Molecule,
    calculator: OpenMMCalculator,
) -> None:
    """Both cages of a built dimer are rigid copies of the monomer."""
    alone = calculator.single_point(cage)
    cage_1, cage_2 = calculator.fragment_energies(dimer)
    assert cage_1 == pytest.approx(alone, abs=0.5)
    assert cage_2 == pytest.approx(alone, abs=0.5)


def test_dimer_energy_terms_add_up(
    cage: stk.BuildingBlock,
    dimer: stk.Molecule,
    calculator: OpenMMCalculator,
) -> None:
    """The derived energies follow their definitions."""
    reference = monomer_energies(cage, calculator)
    energies = dimer_energies(dimer, reference, reference, calculator)

    assert energies.e_min <= energies.e_sp
    assert energies.e_int_sp == pytest.approx(
        energies.e_sp - 2 * reference.as_given
    )
    assert energies.e_bind == pytest.approx(
        energies.e_min - 2 * reference.minimised
    )
    assert energies.e_bind == pytest.approx(
        energies.e_int_min + energies.e_def
    )
    columns = energies.as_dict()
    assert columns["E_bind"] == pytest.approx(energies.e_bind)
    assert len(columns) == 12
