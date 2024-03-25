import numpy as np
import pandas as pd

############### NON CUSTOMIZABLE PART ###################################
def save_network_to_long_rep(G_matrix_path, A_matrix_path, save_path):
    G = np.loadtxt(G_matrix_path)
    A = np.loadtxt(A_matrix_path)

    NR = len(G[0])
    NS = len(G)

    df = pd.DataFrame(columns=['species', 'resource', 'value', 'type'])

    if len(A) != NR:
        print("# rows of A does not match # columns of G")
    if len(A[0])!=NS:
        print("# columns of A does not match # rows of G")

    for i in range(NS):
        for mu in range(NR):
            if A[mu][i]>0:
                df = df.append({
                'species': 'Sp_'+str(i+1),
                'resource': 'Res_'+str(mu+1),
                'value': A[mu][i],
                'type': 'releases'
                }, sort=False, ignore_index=True)
            if G[i][mu] > 0:
                df = df.append({
                'species': 'Sp_'+str(i+1),
                'resource': 'Res_'+str(mu+1),
                'value': G[i][mu],
                'type': 'consumes'
                }, sort=False, ignore_index=True)
    df.to_csv(save_path, index=None, sep='\t', header=None)
    return
############### END OF NON CUSTOMIZABLE PART ################################


G_mat_list="matrix_list/G_matrices.in"
G_folder = "optimal_matrices/consumption/Nr25_Nc25_opt"
A_folder = "optimal_matrices/syntrophy/Nr25_Nc25_opt"
A_suffix = "_optimal_alpha.txt"
save_folder = "optimal_matrices/long_representation/both_modified"

for G_mat_path in np.loadtxt(G_mat_list, dtype="U"):
    A_mat_path = A_folder+G_mat_path[len(G_folder)-4:-4]+A_suffix
    save_path = save_folder+'/G'+G_mat_path[len(G_folder)+9:-3]+'csv'
    save_network_to_long_rep(G_mat_path, A_mat_path, save_path)
