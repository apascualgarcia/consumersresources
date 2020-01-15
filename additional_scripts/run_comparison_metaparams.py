import os
import sys
from datetime import datetime
from separate_file import separate_file, remove_files


config_file = sys.argv[1]
m_list = sys.argv[2]
cores = int(sys.argv[3])
command_name = sys.argv[4]
output_name = sys.argv[5]

with open('./config/' + config_file + '.in') as f:
    mylist = f.read().splitlines()
matrix_list = separate_file(m_list, cores)
seed_number = 0
for i in range(0, cores):
    cmd_core = 'nohup sh -c "'
    for config in mylist:
        out_name = output_name + "_" + config + "_" + m_list + "_" + str(i)
        cmd_core += (command_name + ' ')
        cmd_core += 'config/' + config + '.in'
        cmd_core += " path_to_food_matrix=" + \
            matrix_list[i]
        cmd_core += " path_to_save_file=./data_output/" + output_name + ".out"
        cmd_core += " seed_number=" + str(seed_number)
        cmd_core += "| ts \'[%Y-%m-%d %H:%M:%S]\'\"> ./logs/" + output_name+ ".log"
        cmd_core += " 2>./logs/err" + output_name + '.log'
        seed_number += 1
    #cmd_core = cmd_core[:-3]
    cmd_core += '&'
    os.system(cmd_core)
    now = datetime.now()
    dt_string = now.strftime('%d/%m/%Y %H:%M:%S')
    cmd_write = '[' + dt_string + '] ' + cmd_core
    cmd_write = "echo '" + cmd_write + "' >> logs/commands.log"
    os.system(cmd_write)
print("Launched the commands in the background, please check the appropriate logs to see the progress")
