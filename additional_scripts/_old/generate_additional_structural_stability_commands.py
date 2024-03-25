import os
import sys
from datetime import datetime
from separate_file import separate_file, remove_files
import numpy as np
from os import listdir

LOG_NAME='logs/structural_stability_computations_core'

folder_mode='LRI_connectance_of_no_release_when_eat_Nr50_Nc25'
filename ='commands/struct_stab_computations_modified_perturbation.txt'

command = 'build/compute_critical_Delta'
additional_params='type_of_structural_perturbation=1'


# first get all the config files we have to run
config_folder = 'config/structural_stability'
folders=['common_max_syntrophies', 'maximal_own_syntrophies']
if(len(sys.argv)>1):
    folders=[sys.argv[2]]

files = []
files_wo_folder=[]
config_paths=[config_folder+'/'+a for a in folders]
for c in config_paths:
    for file in listdir(c):
        if folder_mode in file:
            files.append(c+'/'+file)
            files_wo_folder.append(file)
# now we separate files among cores
for i in range(len(files)):
    f = files[i]
    log_name="logs/"+files_wo_folder[i][:-3]+'.log'
    err_log_name="logs/err"+files_wo_folder[i][:-3]+'.log'
    command_core= command +' '+f+ ' ' + additional_params
    command_core+=" | ts \'[%Y-%m-%d %H:%M:%S]\' "
    command_core = command_core+' > '+log_name+' 2>'+err_log_name
    command_core = 'echo "'+command_core+'">> '+filename
    os.system(command_core)
