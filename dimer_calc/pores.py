"""Pores module."""

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

    # hosts = (host1, host2)
    # centroids = (centroid1, centroid2)

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
    # assert np.isclose(pore_distances[0], pore_distance1, atol=1e-6, rtol=0)
    # assert np.isclose(pore_distances[1], pore_distance2, atol=1e-6, rtol=0)


def one_pore(stk_molecule):
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
    # Define calculator object.
    calculator = pm.Inflater(bead_sigma=1.2, centroid=host.get_centroid())

    final_result = calculator.get_inflated_blob(host=host)
    pore = final_result.pore
    pore_distances.append(pore.get_mean_distance_to_com())

    return pore_distances


def check_catenane(one_pore, two_pore, cutoff=0.2):
    # print(one_pore,two_pore)
    if two_pore[0] < (one_pore[0] - cutoff) and two_pore[1] < (
        one_pore[0] - cutoff
    ):
        catenane = True
    else:
        catenane = False
    return catenane


# cage_name='G_17'
# dimer_name="Cage_G_17_63_6_aa"
# cage_dimer=stk.BuildingBlock.init_from_file(f"{dimer_name}.mol")
# cage=stk.BuildingBlock.init_from_file(f"{cage_name}.mol")
#
# cage_pore=one_pore(cage_name,cage)
# dimer_pores=np.min(two_pores(dimer_name,cage_dimer,np.array((0, 0, 0)),np.array((0.41359275, -4.41492435, 12.58085734))))
# print(dimer_pores)
