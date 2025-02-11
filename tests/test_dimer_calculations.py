"""A test script."""

import logging
import pathlib

import stk

from dimer_calculations.axes import by_midpoint, by_smiles, remove_common
from dimer_calculations.dimer_generator import DimerGenerator
from dimer_calculations.pores import get_centroids
from dimer_calculations.utils import remove_aldehyde


def main() -> None:
    """Run script."""
    script_dir = pathlib.Path(__file__).parent.absolute()

    filename = script_dir / ".." / "cages" / "CC1.mol"
    molecule = stk.BuildingBlock.init_from_file(filename)
    molecule = molecule.with_centroid([0, 0, 0])

    list_of_vertices, vertice_size = by_smiles(
        molecule=molecule, smiles_string=remove_aldehyde("NCCN")
    )
    facet_axes, centroid_to_facet = by_midpoint(
        vectors=list_of_vertices,
        vertice_size=vertice_size,
        no_vectors_define_facet=3,
        tolerance=0.1,
    )
    list_of_arenes, centroid_to_arene = by_smiles(
        molecule=molecule,
        smiles_string=remove_aldehyde("O=Cc1cc(C=O)cc(C=O)c1"),
    )
    list_of_windows = remove_common(facet_axes, list_of_arenes)
    centroid_to_window = centroid_to_arene

    dimer_gen = DimerGenerator(
        displacement=10,
        displacement_step_size=1,
        rotation_limit=120,
        rotation_step_size=30,
        slide=False,
    )
    list_of_dimers = dimer_gen.generate(
        molecule=molecule,
        molecule_2=molecule,
        axes_1=list_of_windows[0],
        axes_2=list_of_windows[0],
        displacement_distance=centroid_to_window + centroid_to_window - 2,
    )
    logging.info(len(list_of_dimers))
    logging.info(get_centroids(list_of_dimers[0]["Dimer"]))


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    main()
