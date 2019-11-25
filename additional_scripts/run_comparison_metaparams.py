import os
config_file = 'list_of_config_files_comparison'
matrix_list = 'matrix_list_NR25_NS25'
cores = 6
with open('./config/' + config_file + '.in') as f:
    mylist = f.read().splitlines()

for config in mylist:
    for i in range(1, cores + 1):
        cmd = 'nohup ./build/compute_critical_Delta_matrices '
        cmd += 'config/' + config + '.in'
        cmd += " path_to_food_matrix=./config/" + \
            matrix_list + "_" + str(i) + '.in'
        cmd += " path_to_save_file=./data_output/" + \
            config + "_" + str(i) + ".out"
        cmd += " > ./logs/" + config + "_" + str(i) + ".log"
        cmd += " 2>./logs/err" + config + "_" + str(i) + ".log &"
        os.cmd(cmd)
