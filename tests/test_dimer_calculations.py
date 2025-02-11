"""A test script."""

import pathlib

import stk

from dimer_calculations import Optimiser_functions, utils


def main() -> None:
    """Run script."""
    script_dir = pathlib.Path(__file__).parent.absolute()

    filename = script_dir / ".." / "cages" / "CC1.mol"
    molecule = stk.BuildingBlock.init_from_file(filename)
    molecule = molecule.with_centroid([0, 0, 0])

    list_of_vertices, vertice_size = Optimiser_functions.Axes.BySmiles(
        molecule, smiles_string=utils.remove_aldehyde("NCCN")
    )
    facet_axes, centroid_to_facet = Optimiser_functions.Axes.ByMidpoint(
        molecule,
        vectors=list_of_vertices,
        vertice_size=vertice_size,
        no_vectors_define_facet=3,
        tolerance=0.1,
    )
    list_of_arenes, centroid_to_arene = Optimiser_functions.Axes.BySmiles(
        molecule, smiles_string=utils.remove_aldehyde("O=Cc1cc(C=O)cc(C=O)c1")
    )
    list_of_windows = Optimiser_functions.Axes.RemoveCommon(
        molecule, facet_axes, list_of_arenes
    )
    centroid_to_window = centroid_to_arene

    list_of_dimers = Optimiser_functions.DimerGenerator.generate(
        molecule,
        molecule_2=molecule,
        axes_1=list_of_windows[0],
        axes_2=list_of_windows[0],
        displacement_distance=centroid_to_window + centroid_to_window - 2,
        displacement=10,
        displacement_step_size=1,
        rotation_limit=120,
        rotation_step_size=30,
        slide=False,
    )
    print(len(list_of_dimers))


if __name__ == "__main__":
    main()
