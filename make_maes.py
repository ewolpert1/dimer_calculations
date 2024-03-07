#from Functions import *
#from read_energies_catenated import *
#
#import config
#import utils
import dimer_calc
from dimer_calc.Optimiser_functions import *
from dimer_calc.utils import *
from dimer_calc.cage import *
from dimer_calc.config import *
from dimer_calc.read_energies_catenated_single import *
from dimer_calc.collect_lowest import *
cages_folder_path = config.CAGES_FOLDER_PATH
SCHRODINGER_PATH=config.SCHRODINGER_PATH
GULP_PATH = config.GULP_PATH
XTB_PATH = config.XTB_PATH

#Creates the structures
for filename in os.listdir('cages'):
    if (filename.endswith('.mol')):
        cage_path = os.path.join('cages', filename)
        cage = stk.BuildingBlock.init_from_file(cage_path)
        #gulp_optimizer = GulpDimerOptimizer(gulp_path=GULP_PATH)
        xtb_optimizer = XTBDimerOptimizer(xtb_path=XTB_PATH)
        xtb_optimizer.run_a_cage(filename[:-4], cage, utils.remove_aldehyde('O=Cc1cc(C=O)cc(C=O)c1'), utils.remove_aldehyde('NCCN'))
        #OPLS_optimizer = OPLSDimerOptimizer(SCHRODINGER_PATH=SCHRODINGER_PATH)
        #gulp_optimizer.run_a_cage(filename[:-4], cage, utils.remove_aldehyde('O=Cc1cc(C=O)cc(C=O)c1'), utils.remove_aldehyde('NCCN'))
        #OPLS_optimizer.run_a_cage(filename[:-4], cage, utils.remove_aldehyde('O=Cc1cc(C=O)cc(C=O)c1'), utils.remove_aldehyde('NCCN'))
        
        #run_a_cage(filename[:-4], cage, utils.remove_aldehyde('O=Cc1cc(C=O)cc(C=O)c1'), utils.remove_aldehyde('NCCN'))

        #print('Dimers made now starting to convert files and run calculations...')
#
        #utils.run_a_cage_script(filename[:-4])
#
        #current_dir = os.getcwd()
        #cage_number=f'Cage{filename[:-4]}'
        #minimum_energy = read_energy(f'{cage_number}_mae', cage_number,f"{current_dir}","lowest_dimers")
        #run_calc(current_dir,cage_number,"lowest_dimers")


        