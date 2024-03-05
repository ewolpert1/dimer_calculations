from Functions import *
#from read_energies_catenated import *
from collect_lowest import *
import config
import utils
from read_energies_catenated_single import *

cages_folder_path = config.CAGES_FOLDER_PATH
SCHRODINGER_PATH=config.SCHRODINGER_PATH

#Creates the structures
for filename in os.listdir('cages'):
    if (filename.endswith('.mol')):
        cage_path = os.path.join('cages', filename)
        cage = stk.BuildingBlock.init_from_file(cage_path)
        run_a_cage(filename[:-4], cage, utils.remove_aldehyde('O=Cc1cc(C=O)cc(C=O)c1'), utils.remove_aldehyde('NCCN'))

        print('Dimers made now starting to convert files and run calculations...')

        utils.run_a_cage_script(filename[:-4])

        current_dir = os.getcwd()
        cage_number=f'Cage{filename[:-4]}'
        minimum_energy = read_energy(f'{cage_number}_mae', cage_number,f"{current_dir}","lowest_dimers")
        run_calc(current_dir,cage_number,"lowest_dimers")
