import stk

from dimer_calc.cage import *
from dimer_calc.Optimiser_functions import *


class Dimer:
    molecule: stk.Molecule


filename = "cages/G_17.mol"
molecule = stk.BuildingBlock.init_from_file(filename)

"""
This function returns the vectors between the centroid of the molecule and the smiles string and the magnitude/distance.

You need to give it:
    molecule: stk.Molecule
    smiles_string: str. Here there is a remove_aldehyde function you can use which takes away the aldehydes so you can use the precursor string exactly

"""


list_of_vertices, vertice_size = Axes.BySmiles(
    molecule, smiles_string=utils.remove_aldehyde("NCCN")
)
print(vertice_size)

"""
You can then take these vectors and use them to determine the centroids of the facets, even if theyre windows.

For this you need to give it:
    molecule: stk.Molecule
    vectors: list of vertices
    vertice_size: float
    no_vectors_define_facet: int. This is the number of vectors that define a facet. For example for a 6+4 cage it is 3 but for truncated octahedron it could be 6
    tolerance: float. As this function measures the distance between the vertices to get the ones forming a facet, this is the tolerance for similar the distances have to be where 0.1 is 10% (I think).
"""

facet_axes, facet_size = Axes.ByMidpoint(
    molecule,
    vectors=list_of_vertices,
    vertice_size=vertice_size,
    no_vectors_define_facet=3,
    tolerance=0.1,
)  # if this isnt 4+6 then facet is just windows

list_of_arenes, arene_size = Axes.BySmiles(
    molecule, smiles_string=utils.remove_aldehyde("O=Cc1cc(C=O)cc(C=O)c1")
)

"""
This function is for when you have a symmetry lower than your molecule e.g. 4+6 so you can treats the windows and the arenes separately.

"""

list_of_windows = Axes.RemoveCommon(molecule, facet_axes, list_of_arenes)

window_size = arene_size

"""
Once your axes are defined this function generates the dimers.

For this you need:
    molecule: stk.Molecule
    axes: The axis you want to displace along from the orginal cage
    second_cage_orientation: The orientation of the displaced cage. This funciton rotates the molecule so the vector you give is facing the original cage. For window to arene you just need it the same, but for window-to-window or window-to-arene the cage must be inverted, hence why it is -axis.
    displacement_distance: float. This is the starting distance of the displacement between the two molecules. (It actually starts 2 angstroms closer)
    displacement: float. This+displacement_distance is the maximum distance of displacement
    displacement_step_size: float. This is the step size of the displacement.
    rotation_limit: float. This is the maximum rotation of the displaced cage.
    rotation_step_size: float. This is the step size of the rotation.
    slide: bool. This is whether you want translating ove the face or not.


It returns a dictionary with the following information:

    'Displacement shell': The "shell" of the displacement, i.e. which "displacement" value you're doing
    'Slide': Which position along the face you're at if you're sliding, doesnt actually equate to where but 0 is on axis and not 0 is off axis
    'Rotation': The rotation of the cage
    'Displacement centroid': The actual displacement of the centroid of the neighbouring cage
    'Dimer': the stk.Molecule of the dimer
"""

list_of_dimers = DimerGenerator.generate(
    molecule,
    list_of_windows[0],
    -list_of_windows[0],
    displacement_distance=window_size,
    displacement=5,
    displacement_step_size=1,
    rotation_limit=120,
    rotation_step_size=30,
    slide=False,
)

"""

If you want to constrain the atom positions this function will return the atom positions of the cage you want to fix. It constrains the smiles string unless metal_atom is not none.

For this you need:
    dimer: stk.Molecule
    smiles_string: str. The smiles string of the cage you want to fix
    metal_atom: int. The atom index of the metal atom. If none then it will just fix the smiles string.

I'm pretty sure this only works if your molecule string contains an imine. I imagine you wont use this so wont edit it to be better but lmk if you do want it.

"""

fixed_atom_set = CageOperations.fix_atom_set(
    list_of_dimers[0]["Dimer"], utils.remove_aldehyde("NCCN"), metal_atom=None
)

"""
This is the optimisation part. It will optimise the dimer using GULP, OPLS or XTB.

For this you need:
    dimer_entry: The dictionary entry of the dimer you want to optimise
    output_path: str. Folder you want the optimisation to happen in
    path: path to optimisation software i.e. GULP_PATH, SCHRODINGER_PATH, XTB_PATH
    fixed_atom_set: set. The set of atoms you want to fix, assumed none if not there.

This is what you'll need to add your own CGOptimiser or something to. You can do the in the DimerOptimizer class in the Optimiser_functions.py file.

"""


for dimer_entry in list_of_dimers:
    # gulp_optimizer=DimerOptimizer.optimise_dimer_gulp(dimer_entry['Dimer'],f"gulp_shell_{dimer_entry['Displacement shell']}_slide_{dimer_entry['Slide']}_rot_{dimer_entry['Rotation']}",GULP_PATH,fixed_atom_set)
    # OPLS_optimizer=DimerOptimizer.optimise_dimer_OPLS(dimer_entry['Dimer'],f"OPLS_shell_{dimer_entry['Displacement shell']}_slide_{dimer_entry['Slide']}_rot_{dimer_entry['Rotation']}",SCHRODINGER_PATH,fixed_atom_set)
    XTB_optimizer = DimerOptimizer.optimise_dimer_XTB(
        dimer_entry["Dimer"],
        f"XTB_shell_{dimer_entry['Displacement shell']}_slide_{dimer_entry['Slide']}_rot_{dimer_entry['Rotation']}",
        XTB_PATH,
        fixed_atom_set,
    )
