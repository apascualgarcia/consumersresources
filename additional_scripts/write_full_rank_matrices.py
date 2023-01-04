import numpy as np
from os import listdir
from os.path import isfile, join

matrix_folder="optimal_matrices/consumption/Nr25_Nc25"
matrices=[f for f in listdir(matrix_folder) if isfile(join(matrix_folder, f)) and f.startswith('RandTrix')]

write_file="matrix_list/full_rank_opt_consumption_mat_NR25_NS25.in"
f = open(write_file, 'w+')
for mat in matrices:
    if(mat[0:3]=='Ran'):
        M = np.loadtxt(matrix_folder+'/'+mat)
        if(np.linalg.matrix_rank(M)==len(M)):
            f.write(matrix_folder+'/'+mat+'\n')
f.close()
