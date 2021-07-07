import numpy as np
import pandas as pd

G_matrix_path = "optimal_matrices/consumption/Nr25_Nc25/RandTrix_Nr25_Nc25_Nest0.2_Conn0.168.txt"
A_matrix_path = "optimal_matrices/syntrophy/Nr25_Nc25/RandTrix_Nr25_Nc25_Nest0.2_Conn0.168_optimal_alpha.txt_unconstrained"
save_path = "optimal_matrices/long_representation/G_Nest0.2_Conn0.168.csv"



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
