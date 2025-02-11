import pore_mapper as pm


def two_pores(
    # write_name,
    stk_molecule,
    centroid1,
    centroid2,
):
    # Read in host from xyz file.
    host = pm.Host(
        atoms=(
            pm.Atom(id=i.get_id(), element_string=i.__class__.__name__)
            for i in stk_molecule.get_atoms()
        ),
        position_matrix=stk_molecule.get_position_matrix(),
    )

    # Number of centroids can be any number:
    pore_distances = []
    centroids = (centroid1, centroid2)

    for centroid in centroids:
        # Define calculator object.
        calculator = pm.Inflater(bead_sigma=1.2, centroid=centroid)

        # Run calculator on host object, analysing output.
        final_result = calculator.get_inflated_blob(host=host)
        pore = final_result.pore
        pore_distances.append(pore.get_mean_distance_to_com())

    return pore_distances


def one_pore(
    stk_molecule,
):
    # Read in host from xyz file.
    host = pm.Host(
        atoms=(
            pm.Atom(id=i.get_id(), element_string=i.__class__.__name__)
            for i in stk_molecule.get_atoms()
        ),
        position_matrix=stk_molecule.get_position_matrix(),
    )

    # Number of centroids can be any number:
    pore_distances = []
    calculator = pm.Inflater(bead_sigma=1.2, centroid=host.get_centroid())

    final_result = calculator.get_inflated_blob(host=host)
    pore = final_result.pore
    pore_distances.append(pore.get_mean_distance_to_com())

    return pore_distances


def check_catenane(one_pore, two_pore, cutoff=0.2):
    if two_pore[0] < (one_pore[0] - cutoff) and two_pore[1] < (
        one_pore[0] - cutoff
    ):
        catenane = True
    else:
        catenane = False
    return catenane
