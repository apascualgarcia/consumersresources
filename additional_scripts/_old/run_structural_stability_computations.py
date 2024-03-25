import os
import sys
from datetime import datetime
from separate_file import separate_file, remove_files
import numpy as np
from os import listdir
import random
import string

CORES=int(sys.argv[1])

LOG_NAME='logs/structural_stability_computations_core'
command = 'build/compute_critical_Delta'
additional_params=''

# first get all the config files we have to run
config_folder = 'config/structural_stability/NR50_NS25'
folders=['common_max_syntrophies', 'maximal_own_syntrophies', 'no_syntrophy']
if(len(sys.argv)>2):
    folders=[sys.argv[2]]
    print('Using custom input of folder:', folders)
    if len(sys.argv)> 3:
        for i in range(3, len(sys.argv)):
            additional_params+=sys.argv[i]
            additional_params+=' '

files = []
config_paths=[config_folder+'/'+a for a in folders]
for c in config_paths:
    for file in listdir(c):
        files.append(c+'/'+file)

# now we separate files among cores
files_per_core=np.array_split(files, CORES)
for i in range(CORES):
    rnd=''.join([random.choice(string.ascii_letters + string.digits) for n in range(3)])
    log_name=LOG_NAME+'_'+str(i)+'_'+rnd+'.log'
    command_core = '"'
    for j in range(len(files_per_core[i])):
        command_core+= command +' '+files_per_core[i][j] + ' ' + additional_params
        command_core+=" | ts \'[%Y-%m-%d %H:%M:%S]\' "
        command_core+=' && '
    command_core=command_core[:-4]+'"'
    command_core = "nohup sh -c "+command_core+' > '+log_name+' 2>&1 &'
    os.system(command_core)
print("Launched "+str(len(files))+" runs on "+str(CORES)+" cores. Please check appropriate log folders if needed.")
