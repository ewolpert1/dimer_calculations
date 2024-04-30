
import stk
from dimer_calc.Optimiser_functions import *
from dimer_calc.cage import *



#cages_folder_path = config.CAGES_FOLDER_PATH
#SCHRODINGER_PATH=config.SCHRODINGER_PATH
#GULP_PATH = config.GULP_PATH
#XTB_PATH = config.XTB_PATH

class Dimer:
    molecule: stk.Molecule


#def main():
filename='cages/G_17.mol'
#filename_xyz='cages/G_17.xyz'
molecule=stk.BuildingBlock.init_from_file(filename)

'''
Currently BySmiles is the only one that works
'''


list_of_vertices,vertice_size = Axes.BySmiles(molecule,smiles_string=utils.remove_aldehyde('NCCN'))

print(vertice_size)

facet_axes,facet_size =Axes.ByMidpoint(molecule,
    vectors=list_of_vertices,
    vertice_size=vertice_size,
    no_vectors_define_facet=3,
    tolerance=0.1)#if this isnt 4+6 then facet is just windows

list_of_arenes,arene_size = Axes.BySmiles(molecule,smiles_string=utils.remove_aldehyde('O=Cc1cc(C=O)cc(C=O)c1'))

list_of_windows =Axes.RemoveCommon(molecule,facet_axes, list_of_arenes)

window_size=arene_size

list_of_dimers = DimerGenerator.generate(molecule,
        list_of_windows[0],
        -list_of_windows[0],
        displacement_distance=window_size,
        displacement=5,
        displacement_step_size=1,
        rotation_limit=120,
        rotation_step_size=30,
        slide=False,
        )

#print(list_of_dimers[1])
#print(list_of_dimers[0][0]['Dimer'])

fixed_atom_set = CageOperations.fix_atom_set(list_of_dimers[0]['Dimer'],utils.remove_aldehyde('NCCN'), metal_atom=None)

for dimer_entry in list_of_dimers:
    print(dimer_entry['Dimer'])

    #gulp_optimizer = GulpDimerOptimizer(gulp_path=GULP_PATH)
#    gulp_optimizer=DimerOptimizer.optimise_dimer_gulp(dimer_entry['Dimer'],f"gulp_shell_{dimer_entry['Displacement shell']}_slide_{dimer_entry['Slide']}_rot_{dimer_entry['Rotation']}",GULP_PATH,fixed_atom_set)
    #gulp_optimizer=DimerOptimizer.optimise_dimer_gulp(dimer_entry['Dimer'],f"gulp_shell_{dimer_entry['Displacement shell']}_slide_{dimer_entry['Slide']}_rot_{dimer_entry['Rotation']}",GULP_PATH,fixed_atom_set)
    #OPLS_optimizer=DimerOptimizer.optimise_dimer_OPLS(dimer_entry['Dimer'],f"OPLS_shell_{dimer_entry['Displacement shell']}_slide_{dimer_entry['Slide']}_rot_{dimer_entry['Rotation']}",SCHRODINGER_PATH,fixed_atom_set)
    XTB_optimizer=DimerOptimizer.optimise_dimer_XTB(dimer_entry['Dimer'],f"XTB_shell_{dimer_entry['Displacement shell']}_slide_{dimer_entry['Slide']}_rot_{dimer_entry['Rotation']}",XTB_PATH,fixed_atom_set)





#
#    for dimer in list_o_dimers:
#        # Optimise.
#        opt_molecule = stko.Somehting/dimers.optimiser(settings).optimize(dimer.molecule)
#        # Energy.
#        get the energy somehow.
#

    # Table of energy vs. dimer-position?


#if __name__ == '__main__':
#    main()
#
#
##Creates the structures
#for filename in os.listdir('cages'):
#    if (filename.endswith('G_17.mol')):
#        cage_path = os.path.join('cages', filename)
#        cage = stk.BuildingBlock.init_from_file(cage_path)
#        gulp_optimizer = GulpDimerOptimizer(gulp_path=GULP_PATH)
#        #xtb_optimizer = XTBDimerOptimizer(xtb_path=XTB_PATH)
#        #xtb_optimizer.run_a_cage(filename[:-4], cage, utils.remove_aldehyde('O=Cc1cc(C=O)cc(C=O)c1'), utils.remove_aldehyde('NCCN'))
#        #OPLS_optimizer = OPLSDimerOptimizer(SCHRODINGER_PATH=SCHRODINGER_PATH)
#        gulp_optimizer.run_a_cage(filename[:-4], cage, utils.remove_aldehyde('O=Cc1cc(C=O)cc(C=O)c1'), utils.remove_aldehyde('NCCN'))
#        #OPLS_optimizer.run_a_cage(filename[:-4], cage, utils.remove_aldehyde('O=Cc1cc(C=O)cc(C=O)c1'), utils.remove_aldehyde('NCCN'))
#
#        #run_a_cage(filename[:-4], cage, utils.remove_aldehyde('O=Cc1cc(C=O)cc(C=O)c1'), utils.remove_aldehyde('NCCN'))
#
#        #print('Dimers made now starting to convert files and run calculations...')
##
#        #utils.run_a_cage_script(filename[:-4])
##
#        #current_dir = os.getcwd()
#        #cage_number=f'Cage{filename[:-4]}'
#        #energy_reader = XTBEnergyReader(cage_name_xtb=f'{cage_number}_xtb',
#        #                         cage_name=cage_number,
#        #                         current_directory=f"{current_dir}",
#        #                         destination_folder_end="lowest_dimers_xtb")
#        #minimum_energy = energy_reader.read_energy()
#
#        #run_calc(current_dir,cage_number,"lowest_dimers")


