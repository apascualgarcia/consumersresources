import os
import numpy as np


def separate_file(conf_file, CORES):
    with open('matrix_list/' + conf_file + '.in') as f:
        mylist = f.read().splitlines()
    indices_tables = np.array_split(np.arange(0, len(mylist)), CORES)
    filenames = []
    for i in range(len(indices_tables)):
        filenames.append('./config/' + conf_file + '_' + str(i) + '.in')
        f = open(filenames[i], 'w+')
        for j in indices_tables[i]:
            f.write(str(mylist[j]) + '\n')
        f.close()
    return filenames


def remove_files(mylist):
    for file in mylist:
        cmd = "rm " + file
        os.system(cmd)
