import os
import sys
from datetime import datetime
from separate_file import separate_file, remove_files

config_file = sys.argv[1]
matrix_list = sys.argv[2]
cores = int(sys.argv[3])

with open('./config/' + config_file + '.in') as f:
    mylist = f.read().splitlines()
matrix_list = separate_file(matrix_list, cores)
seed_number = 0
for i in range(0, cores):
    cmd_core = 'nohup '
    for config in mylist:
        cmd_core += './build/compute_critical_Delta_matrices '
        cmd_core += 'config/' + config + '.in'
        cmd_core += " path_to_food_matrix=./config/" + \
            matrix_list[i] + '.in'
        cmd_core += " path_to_save_file=./data_output/" + \
            config + "_" + str(i) + ".out"
        cmd_core += " seed_number=" + str(seed_number)
        cmd_core += " | ts '[%Y-%m-%d %H:%M:%S]' > ./logs/" + \
            config + "_" + str(i) + ".log"
        cmd_core += " 2>./logs/err" + config + "_" + \
            str(i) + ".log && "
        seed_number += 1
    cmd_core = cmd_core[:-3]
    cmd_core = cmd_core + '&'
    print(cmd_core)
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    cmd_write = '[' + dt_string + '] ' + cmd_core
    cmd_write = 'echo "' + cmd_write + '" >> logs/commands.log'
    os.system(cmd_write)
print("Launched the commands in the background, please check the appropriate logs to see the progress")
