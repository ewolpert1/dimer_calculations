# Introduction
Hi Harry and Andrew, welcome to my dimer calculations code. I thought it would be best to get this out to you ASAP. Due to this and my distinct lack of coding skills this is very much in beta mode and has not been thoroughly tested. Therefore I'm sure some issues will arise and I imagine the usability isn't great and there are definitely dead functions, sorry about that, but any thoughts/suggestions of what could help would be much appreciated!

The code works by taking cage molecules and creates dimers by displacing a cage along each of the high symmetry axes, which are: window to window (ww), window to arene (wa), window to vertex (wa), and arene to arene (aa). It displaces each of the cages such that the facets defined are 1 angstrom away (e.g. the window to window facet) and then displaces it at 1 angstrom increments up to 8 angstroms away. It also displaces each of the cages around the facet to properly sample configurational space, and rotates them about the high symmetry axis in 30 degree increments. The values of each of these can be changed. Currently its hard coded in as it would need to be done within the optimiser_functions.py file but watch this space for more functionality in the future. The code checks for overlap between the cages, as well as catenation and removes any cages with overlap and/or catenation. These are set with tolerance factors which again you can change but are currently hard coded into it. 

# Optimisers

There are three different optimisers you have the pleasure of choosing between:

OPLS, GULP, and XTB whose paths need to be defined in the main script (or the config.py file). I have included an example of how to run this code in the make_dimers.py file. Each of the optimisers must constain the following variables:

1. The name of the cage
2. The cage itself (as an stk.BuildingBlock)
3. The trialdehyde smile string (the string defining the "arene" facet)
4. The diamine smile string (the string defining the vertices of the cage)
5. The symmetry of the cage. Currently there are three options but I'm updating these as and when theyre useful. The options are:
    sym="4_6" which is the 4+6 topology with tetrahedral symmetry (classic cage 3)
    sym="oct" octahedral cage but with octahedral symmetry (like WaaF-1)
    sym="trunc_oct" truncated octahedral symmetry (Like Mastarlerzs cage)

There are also other variables such as metal atoms, paralleisations, choice in optimisation or not, choice in constraints or not. But these havent been fully implemented yet so I would avoid using for now. If you do want to use a metal atom that is not Rh, let me know and I will bump that up my to do list.

Constraints mean that I am constraining the atoms which make up the vertices of the cages. I think this only works when when you have a diamine or a metal atom in your system as it either constrains the atoms neighbouring the imine or the metal atom itself. This is so you dont get much rearrangement in the optimisation and also to conserve the shape of the cage.

With that so far you should be able to create all the scripts necessary to run the dimer calculations for OPLS, XTB, and GULP. If you cant please tell me and I will update the the README/scripts with better instructions/functionality.

# Running scripts

You will see when running the script it creates two folders: Cage{name} and Cage{name}_XXX where XXX differs depending on which optimiser you are using. The first folder contains all the .mol files of the dimers you create into four folders based on the type of packing. They're named as "Cage_{name}_{displacement perpendicular to high sym axis}_{rotation}_{displacement along high sym axis}_aa.mol" where the numbers in it dont correspond to distance/rotation angle. You can have a look at these to make sure its doing what you expect, if not it might be an issue with what smiles string you've used (for example if the backbone of the trialdehyde also occurs in the vertex/diamine). The other folder is where your dimer calculations will be run.

For GULP, the dimer calculations are actually performed as you run the script. 

For OPLS you will have to go into the _mae folder and run the "run_a_cage.sh" script. The thing thats most time consuming about this is the converting mol files to mae files which is done through schrodinger. I imagine there might be a better way of doing this which we should look into, especially for you Harry if you need to do this in a high throughput way. 

For XTB its similar to OPLS where you have to go into the _XTB folder and run "constrain.sh". The thing that is time consuming about this is that its XTB :)

I have also written energy_readers which then read the energy of your outputted structures and collates the lowest energy dimers (within 50kJ/mol of the lowest configuration for each packing type if I'm not mistaken). If the first calculations you did is constrained you can then run a second, unconstrained calculation on the lowest energy dimers. I have energy_readers for GULP and OPLS but not for XTB yet. TBH I havent looked/tested these in a while so I will play around and update this README when I have.

# HARRY

We will probably need to edit this quite a bit. I think if we make a new function i.e. run_a_molecule instead of run_a_cage and then also a new function similar to generate_new_cages but based on the calculated shape of your molecule instead of smiles it could work quite well. I suggest you make a new branch(?) or whatever its called and then add in some functions that will work on your molecules.

Happy dimer calculating!
