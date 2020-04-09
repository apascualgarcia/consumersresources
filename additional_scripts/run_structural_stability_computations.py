import os
import sys
from datetime import datetime
from separate_file import separate_file, remove_files
import numpy as np
from os import listdir

CORES=int(sys.argv[1])
print(CORES)
LOG_NAME='logs/structural_stability_computations_core'

command = 'build/compute_critical_Delta'

# first get all the config files we have to run
config_folder = 'config/structural_stability'
folders=['common_max_syntrophies', 'maximal_own_syntrophies', 'no_syntrophy']

files = []
config_paths=[config_folder+'/'+a for a in folders]
for c in config_paths:
    for file in listdir(c):
        files.append(c+'/'+file)

# now we separate files among cores
files_per_core=np.array_split(files, CORES)
for i in range(CORES):
    log_name=LOG_NAME+'_'+str(i)+'.log'
    command_core = '"'
    for j in range(len(files_per_core[i])):
        command_core+= command +' '+files_per_core[i][j]+' && '
    command_core=command_core[:-4]+'"'
    command_core = "nohup sh -c "+command_core+' > '+log_name+' &'
    os.system(command_core)
