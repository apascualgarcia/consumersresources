import numpy as np
import glob

suffix = "_opt_mat_alpha0=1_intra_specific_syntrophy=allowed_seed_number=*_gamma0=1"
matrix_list = "matrix_list/full_rank_opt_consumption_mat_NR25_NS25.in"
[NR, NS]=[25, 25]

for G_mat in np.loadtxt(matrix_list, dtype='U'):
    match_name = G_mat[:-4]+"_optimal_alpha.txt"+suffix
    mat_files = glob.glob(match_name)
    A_mat = np.zeros((NR, NS))
    number_of_matrices=0
    for mat_path in mat_files:
        to_add = np.loadtxt(mat_path)
        if to_add.size>0:
            A_mat+=to_add
            number_of_matrices+=1
    A_mat/=number_of_matrices
    np.savetxt(G_mat[:-4]+"_optimal_alpha_meta_matrix.txt", A_mat, dtype=np.float64)
