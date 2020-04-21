import os
import sys
from datetime import datetime
from separate_file import separate_file, remove_files
import numpy as np
import os

file_name = sys.argv[1]
CORES = int(sys.argv[2])
commands = np.array([" ".join(f) for f in np.loadtxt(file_name, dtype='U')])

commands_per_core = np.array_split(commands, CORES)

for i in range(CORES):
    core_command = ''
    for j in range(len(commands_per_core[i])):
        core_command+="nohup "+commands_per_core[i][j]+' && '
    core_command=core_command[:-4]+' &'
    print(core_command)
    #os.system(core_command)
