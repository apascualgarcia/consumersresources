import numpy as np

file_list = 'matrix_list/Nr25_Nc25/all_G_matrices.in'
mat_list = np.loadtxt(file_list, dtype='U')

return_string=[]
for mat in mat_list:
    M = np.loadtxt(mat)
    max_rank = min(len(M), len(M[0]))
    rank = np.linalg.matrix_rank(M)
    if rank!=max_rank:
        a=mat.split("_",7)
        nest=a[5][4:]
        conn=a[6][4:-4]
        return_string.append("{"+nest+','+conn+"}")
print(return_string)
