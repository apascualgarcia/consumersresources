import numpy as np
from os import listdir
from os.path import isfile, join

matrix_folder="test_matrices"
matrices=[f for f in listdir(matrix_folder) if isfile(join(matrix_folder, f))]

write_file="matrix_list/new_matrices.in"
f = open(write_file, 'w+')
for mat in matrices:
    if(mat[0:3]=='Ran'):
        M = np.loadtxt(matrix_folder+'/'+mat)
        if(np.linalg.matrix_rank(M)==len(M)):
            f.write(matrix_folder+'/'+mat+'\n')
f.close()
