"""TEst CG model."""

import pathlib

import cgexplore
import matplotlib as mpl
import matplotlib.pyplot as plt
import openmm
import stk

from dimer_calc.Optimiser_functions import Axes, DimerGenerator


def get_openmm_energy(
    molecule: stk.Molecule,
    forcefield: cgexplore.forcefields.ForceField,
    name: str,
    output_directory: pathlib.Path,
) -> float:
    """Compute OpenMM energy."""
    assigned_system = forcefield.assign_terms(
        molecule=molecule,
        name=name,
        output_dir=output_directory,
    )

    # Get energy.
    return (
        cgexplore.optimisation.CGOMMOptimizer(
            fileprefix=name,
            output_dir=output_directory,
            max_iterations=0,
            platform=None,
        )
        .calculate_energy(assigned_system)
        .value_in_unit(openmm.unit.kilojoules_per_mole)
    )


def main() -> None:
    raise SystemExit("delete me")
    working_dir = pathlib.Path(
        "/home/atarzia/workingspace/pore_topology_optimisation/"
    )
    output_dir = working_dir / "dimer_test"
    system_dir = working_dir / "systems"

    systems = (
        "sopt_00051616101830",
        "sopt_1104335528716",
        "sopt_20142345301035",
    )

    # Define the forcefield for all.
    print("remove this when you move to pore topo")
    mbead = cgexplore.molecular.CgBead(
        element_string="Pd",
        bead_class="m",
        bead_type="m",
        coordination=4,
    )
    nbead = cgexplore.molecular.CgBead(
        element_string="C",
        bead_class="n",
        bead_type="n",
        coordination=3,
    )
    bbead = cgexplore.molecular.CgBead(
        element_string="Pb",
        bead_class="b",
        bead_type="b",
        coordination=2,
    )
    abead = cgexplore.molecular.CgBead(
        element_string="Ba",
        bead_class="a",
        bead_type="a",
        coordination=2,
    )
    cbead = cgexplore.molecular.CgBead(
        element_string="Ag",
        bead_class="c",
        bead_type="c",
        coordination=2,
    )

    present_beads = (mbead, nbead, cbead, bbead, abead)
    definer_dict = {
        # Nonbondeds.
        "m": ("nb", 10.0, 1.0),
        "n": ("nb", 10.0, 1.0),
        "a": ("nb", 10.0, 1.0),
        "b": ("nb", 10.0, 1.0),
        "c": ("nb", 10.0, 1.0),
    }
    nonbonded_terms = []

    for interaction_key in definer_dict:
        interaction_list = definer_dict[interaction_key]
        nonbonded_terms.append(
            cgexplore.terms.TargetNonbonded(
                bead_class=interaction_key[0],
                bead_element=next(
                    i.element_string
                    for i in present_beads
                    if i.bead_type == interaction_key[0]
                ),
                epsilon=openmm.unit.Quantity(
                    value=interaction_list[1],
                    unit=openmm.unit.kilojoules_per_mole,
                ),
                sigma=openmm.unit.Quantity(
                    value=interaction_list[2],
                    unit=openmm.unit.angstrom,
                ),
                force="custom-lj",
            )
        )
    forcefield = cgexplore.forcefields.ForceField(
        identifier="test",
        prefix="test",
        present_beads=present_beads,
        bond_targets=(),
        angle_targets=(),
        torsion_targets=(),
        nonbonded_targets=tuple(nonbonded_terms),
        vdw_bond_cutoff=2,
        verbose=False,
    )

    for system, tstr in zip(systems, ("6P12", "8P12", "6P8"), strict=False):
        filename = working_dir / system_dir / f"{system}_optc.mol"
        print(f"doing {system} with {tstr}")
        molecule = stk.BuildingBlock.init_from_file(filename)

        # Determine topology dependant vectors of interest.
        if tstr == "6P12":
            ag_vectors, ag_size = Axes.BySmiles(
                molecule,
                smiles_string="[Ag]",
            )
            hexagon_vectors, hexagon_mean_distance = (
                Axes.ByMidpoint_mindistance(
                    molecule,
                    vectors=ag_vectors,
                    vertice_size=ag_size,
                    no_vectors_define_facet=3,
                    tolerance=0.1,
                )
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
            hexagon_vectors, hexagon_mean_distance = (
                Axes.Square_ByMidpoint_mindistance(
                    molecule,
                    vectors=ag_vectors,
                    vertice_size=ag_size,
                    no_vectors_define_facet=4,
                    tolerance=0.1,
                )
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
                "tab:blue",
            ),
            "hs": (
                hexagon_vectors[0],
                -square_vectors[0],
                hexagon_mean_distance,
                "hex-squ",
                "tab:orange",
            ),
            "ss": (
                square_vectors[0],
                -square_vectors[0],
                square_mean_distance,
                "squ-squ",
                "tab:green",
            ),
        }

        dimers = {}
        datas = {ac: [] for ac in axis_combinations}
        for ac in axis_combinations:
            dimers[ac] = DimerGenerator.generate(
                molecule,
                axes=axis_combinations[ac][0],
                second_cage_orientation=axis_combinations[ac][1],
                displacement_distance=axis_combinations[ac][2],
                displacement=3,
                displacement_step_size=0.5,
                rotation_limit=120,
                rotation_step_size=5,
                slide=False,
                overlap_tolerance=0.95,
            )

            print(f"made {len(dimers[ac])} {ac} dimers")
            print("need to make this use a saving DB")
            for dimer_entry in dimers[ac]:
                dimer_molecule = dimer_entry["Dimer"]

                dimer_name = (
                    f"{system}_{ac}_s_{dimer_entry['Displacement shell']}_"
                    f"sl_{dimer_entry['Slide']}_r_{dimer_entry['Rotation']}"
                )
                dimer_molecule.write(output_dir / f"{dimer_name}.mol")

                # Calculate energy.
                datas[ac].append(
                    (
                        dimer_entry["Displacement shell"],
                        dimer_entry["Rotation"],
                        get_openmm_energy(
                            molecule=dimer_molecule,
                            name=dimer_name,
                            forcefield=forcefield,
                            output_directory=output_dir,
                        ),
                    )
                )

        # Plot energy as a function of
        fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(16, 5))
        vmin = -20
        vmax = 0
        for ac, ax in zip(datas, axs, strict=False):
            print(
                min([i[2] for i in datas[ac]]), max([i[2] for i in datas[ac]])
            )
            ax.scatter(
                [i[0] for i in datas[ac]],
                [i[1] for i in datas[ac]],
                c=[i[2] for i in datas[ac]],
                vmin=vmin,
                vmax=vmax,
                alpha=1.0,
                edgecolor="k",
                s=200,
                marker="s",
                cmap="Blues_r",
            )

            ax.tick_params(axis="both", which="major", labelsize=16)
            ax.set_xlabel("displacement", fontsize=16)
            ax.set_ylabel("rotation", fontsize=16)
            ax.set_title(axis_combinations[ac][3], fontsize=16)

        cbar_ax = fig.add_axes([1.01, 0.2, 0.02, 0.7])
        cmap = mpl.cm.Blues_r
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        cbar = fig.colorbar(
            mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
            cax=cbar_ax,
            orientation="vertical",
        )
        cbar.ax.tick_params(labelsize=16)
        cbar.set_label("$E$", fontsize=16)

        fig.tight_layout()
        fig.savefig(
            output_dir / f"{system}_energies.png",
            dpi=360,
            bbox_inches="tight",
        )
        plt.close()


if __name__ == "__main__":
    main()
