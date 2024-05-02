import stk
from dimer_calc.Optimiser_functions import Axes, DimerGenerator
import pathlib


def main() -> None:
    working_dir = pathlib.Path("/home/atarzia/workingspace/pore_topology_optimisation/")
    output_dir = working_dir / "dimer_test"
    system_dir = working_dir / "systems"

    systems = (
        "sopt_00051616101830",
        "sopt_1104335528716",
        "sopt_20142345301035",
    )

    for system, tstr in zip(systems, ("6P12", "8P12", "6P8")):
        filename = working_dir / system_dir / f"{system}_optc.mol"
        print(f"doing {system} with {tstr}")
        molecule = stk.BuildingBlock.init_from_file(filename)

        # Determine topology dependant vectors of interest.
        if tstr == "6P12":
            ag_vectors, ag_size = Axes.BySmiles(
                molecule,
                smiles_string="[Ag]",
            )
            hexagon_vectors, hexagon_mean_distance = Axes.ByMidpoint_mindistance(
                molecule,
                vectors=ag_vectors,
                vertice_size=ag_size,
                no_vectors_define_facet=3,
                tolerance=0.1,
            )
            square_vectors, square_mean_distance = Axes.BySmiles(
                molecule,
                smiles_string="[Pd]",
            )
        elif tstr == "8P12":
            ag_vectors, ag_size = Axes.BySmiles(
                molecule,
                smiles_string="[Ag]",
            )
            hexagon_vectors, hexagon_mean_distance = Axes.Square_ByMidpoint_mindistance(
                molecule,
                vectors=ag_vectors,
                vertice_size=ag_size,
                no_vectors_define_facet=4,
                tolerance=0.1,
            )

            square_vectors, square_mean_distance = Axes.BySmiles(
                molecule,
                smiles_string="C",
            )

        elif tstr == "6P8":
            continue

        # Visualise the vectors.
        vector_file = output_dir / f"{system}_vects.xyz"
        with vector_file.open("w") as f:
            f.write(f"{len(hexagon_vectors)+len(square_vectors)}\n\n")
            for i in hexagon_vectors:
                f.write(f"C {i[0]} {i[1]} {i[2]}\n")
            for i in square_vectors:
                f.write(f"S {i[0]} {i[1]} {i[2]}\n")

        axis_combinations = {
            "hh": (
                hexagon_vectors[0],
                -hexagon_vectors[0],
                hexagon_mean_distance,
                "hex-hex",
            ),
            "hs": (
                hexagon_vectors[0],
                -square_vectors[0],
                hexagon_mean_distance,
                "hex-squ",
            ),
            "ss": (
                square_vectors[0],
                -square_vectors[0],
                square_mean_distance,
                "squ-squ",
            ),
        }

        dimers = {}
        for ac in axis_combinations:
            dimers[ac] = DimerGenerator.generate(
                molecule,
                axes=axis_combinations[ac][0],
                second_cage_orientation=axis_combinations[ac][1],
                displacement_distance=axis_combinations[ac][2],
                displacement=3,
                displacement_step_size=0.5,
                rotation_limit=120,
                rotation_step_size=15,
                slide=False,
                overlap_tolerance=0.95,
            )

            print(f"made {len(dimers[ac])} {ac} dimers")

            for dimer_entry in dimers[ac]:
                dimer_molecule = dimer_entry["Dimer"]
                dimer_unopt_file = (
                    f"{system}_{ac}_s_{dimer_entry['Displacement shell']}_"
                    f"sl_{dimer_entry['Slide']}_r_{dimer_entry['Rotation']}"
                )
                dimer_molecule.write(output_dir / f"{dimer_unopt_file}.mol")


if __name__ == "__main__":
    main()
