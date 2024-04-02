import numpy as np
from os import listdir
from os.path import isfile, join
import re



matrix_root_folder="./matrices/Nr25_Nc25/syntrophy/optimized_alpha=0.005_gamma=0.5/Binary"
seeds = range(50)
for seed in seeds:
    matrix_folder = matrix_root_folder+"/seed_"+str(seed)
    matrices=[matrix_folder+'/'+f for f in listdir(matrix_folder) if isfile(join(matrix_folder, f)) and not f.startswith('.')]
    for mat in matrices:
        start_ending = mat.find(".txt")
        new_name = mat[:start_ending+4]
        if mat[-6:] == "energy":
            new_name = new_name[:-4]+"_energy"
        f = open(mat)
        header_1 = f.readline()
        header_2 = f.readline()
        np.savetxt(new_name, np.loadtxt(mat), fmt='%d', header=header_1 + " "+ header_2)
