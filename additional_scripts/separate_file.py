import os
import numpy as np

config_folder = './config'
config_file = 'matrix_list_NR25_NS25'
CORES = 20


def separate_file(conf_file):
    with open(conf_file) as f:
        mylist = f.read().splitlines()
    indices_tables = np.array_split(np.arange(0, len(mylist)), CORES)
    filenames = []
    for i in range(len(indices_tables)):
        filenames.append('./config/' + config_file + '_' + str(i) + '.in')
        f = open(filenames[i], 'w+')
        for j in indices_tables[i]:
            f.write(str(mylist[j]) + '\n')
        f.close()
    return filenames


def remove_files(mylist):
    for file in mylist:
        cmd = "rm " + file
        os.system(cmd)


config_files = separate_file(config_folder + '/' + config_file + '.in')
remove_files(config_files)
