import os
import sys
from datetime import datetime
from separate_file import separate_file, remove_files
import numpy as np
import os
import random

file_name = sys.argv[1]
CORES = int(sys.argv[2])
command_list = np.loadtxt(file_name, dtype='U')
if command_list.ndim==1:
    commands = [" ".join(command_list)]
else:
    commands = np.array([" ".join(f) for f in command_list])
np.random.shuffle(commands)
commands_per_core = np.array_split(commands, CORES)

log_name=str(file_name).split("/")
log_name=log_name[1].split(".")[0]


for i in range(CORES):
    core_command = '('
    for j in range(len(commands_per_core[i])):
        log_file = "logs/"+log_name+"_core_"+str(i)+".log"
        err_file = "logs/err_"+log_name+"_core_"+str(i)+".log"

        core_command+='nohup sh -c "'+commands_per_core[i][j]
        core_command+="| ts '[%Y-%m-%d %H:%M:%S]' "
        core_command+='">'+log_file+" 2>"+err_file
        core_command+="; wait;"
    core_command=core_command[:-7]+') &'
    print(core_command)
    os.system(core_command)
