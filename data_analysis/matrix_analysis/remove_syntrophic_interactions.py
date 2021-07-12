import numpy as np



def remove_iss_links(G_matrix_path, A_matrix_path):
    G = np.loadtxt(G_matrix_path)
    A = np.loadtxt(A_matrix_path)

    for i in range(len(G)):
        for mu in range(len(G[i])):
            if G[i][mu]*A[mu][i] > 0:
                A[mu][i]=0
    np.savetxt(A_matrix_path[:-4]+"_no_iss.txt", A, fmt='%i')
    return


G_matrix_path = "optimal_matrices/consumption/Nr25_Nc25/RandTrix_Nr25_Nc25_Nest0.4_Conn0.1232.txt"
A_matrix_path = "optimal_matrices/consumption/Nr25_Nc25/RandTrix_Nr25_Nc25_Nest0.4_Conn0.1232_optimal_alpha.txt"

G_mat_list="matrix_list/full_rank_opt_consumption_mat_NR25_NS25.in"
G_folder = "optimal_matrices/consumption/Nr25_Nc25"
A_folder = "optimal_matrices/syntrophy/Nr25_Nc25/9Jul21"
A_suffix = "_optimal_alpha.txt_optimize_matrices_alpha0=1_intra_specific_syntrophy=allowed_gamma0=1"

for G_mat_path in np.loadtxt(G_mat_list, dtype="U"):
    A_mat_path = A_folder+G_mat_path[len(G_folder):-4]+A_suffix
    remove_iss_links(G_mat_path, A_mat_path)
