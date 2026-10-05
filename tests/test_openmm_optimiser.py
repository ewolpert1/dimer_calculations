"""Tests for the OpenMM optimiser. Skipped when OpenMM or OpenFF is missing."""

import pathlib

import numpy as np
import pytest
import stk

pytest.importorskip("openmm")
pytest.importorskip("openff.toolkit")
pytest.importorskip("openff.interchange")

from dimer_calculations.axes import (
    by_midpoint,
    by_smiles,
    remove_common,
)
from dimer_calculations.cage import Cage
from dimer_calculations.dimer_generator import DimerGenerator
from dimer_calculations.optimiser_functions import (
    OpenMMDimer,
    optimise_dimer_openmm,
)
from dimer_calculations.utils import remove_aldehyde

CC1 = pathlib.Path(__file__).parent.parent / "cages" / "CC1.mol"


@pytest.fixture(scope="module")
def dimer() -> stk.Molecule:
    """Build one window-to-window CC1 dimer, as in the smoke test."""
    cage = stk.BuildingBlock.init_from_file(CC1).with_centroid([0, 0, 0])
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


def test_optimise_dimer_openmm_writes_the_optimised_dimer(
    dimer: stk.Molecule, tmp_path: pathlib.Path
) -> None:
    """The optimised dimer and its energies are written."""
    output_dir = tmp_path / "OpenMM_shell_0"
    optimise_dimer_openmm(
        dimer=dimer, output_dir=str(output_dir), max_iterations=200
    )

    optimised = stk.BuildingBlock.init_from_file(
        tmp_path / "OpenMM_shell_0_opt.mol"
    )
    assert optimised.get_num_atoms() == dimer.get_num_atoms()
    assert not np.allclose(
        optimised.get_position_matrix(), dimer.get_position_matrix()
    )

    lines = (output_dir / "openmm_opt.out").read_text().splitlines()
    values = dict(line.split(": ") for line in lines)
    assert float(values["final_energy_kJ_mol"]) < float(
        values["initial_energy_kJ_mol"]
    )


def test_fixed_atoms_keep_their_positions(dimer: stk.Molecule) -> None:
    """Fixed atoms do not move; the rest of the dimer relaxes."""
    fixed_atom_set = Cage.fix_atom_set(dimer, "NCCN")
    assert fixed_atom_set

    optimised = OpenMMDimer(max_iterations=200).optimize(
        mol=dimer, fixed_atom_set=fixed_atom_set
    )
    shift = np.linalg.norm(
        optimised.get_position_matrix() - dimer.get_position_matrix(),
        axis=1,
    )
    free = np.ones(len(shift), dtype=bool)
    free[fixed_atom_set] = False
    assert shift[fixed_atom_set].max() < 1e-6
    assert shift[free].max() > 0.01


def test_a_finished_dimer_is_skipped(
    dimer: stk.Molecule, tmp_path: pathlib.Path
) -> None:
    """An existing ``_opt.mol`` is left untouched."""
    output_dir = tmp_path / "OpenMM_shell_0"
    done = tmp_path / "OpenMM_shell_0_opt.mol"
    done.write_text("finished earlier")
    optimise_dimer_openmm(dimer=dimer, output_dir=str(output_dir))
    assert done.read_text() == "finished earlier"


def test_fixed_atom_ids_outside_the_dimer_are_rejected(
    dimer: stk.Molecule,
) -> None:
    """A wrong atom id raises an error instead of being ignored."""
    with pytest.raises(ValueError, match="outside"):
        OpenMMDimer().optimize(
            mol=dimer, fixed_atom_set=[dimer.get_num_atoms()]
        )
