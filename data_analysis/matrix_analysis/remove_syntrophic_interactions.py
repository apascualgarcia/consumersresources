import numpy as np

G_matrix_path = "optimal_matrices/consumption/Nr25_Nc25/RandTrix_Nr25_Nc25_Nest0.4_Conn0.1232.txt"
A_matrix_path = "optimal_matrices/consumption/Nr25_Nc25/RandTrix_Nr25_Nc25_Nest0.4_Conn0.1232_optimal_alpha.txt"

G = np.loadtxt(G_matrix_path)
A = np.loadtxt(A_matrix_path)

for i in range(len(G)):
    for mu in range(len(G[i])):
        if G[i][mu]*A[mu][i] > 0:
            A[mu][i]=0
np.savetxt(A_matrix_path[:-4]+"_no_iss.txt", A, fmt='%i')
