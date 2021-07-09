import numpy as np
import pandas as pd
import consumer_resource_data_analysis as cf

G_matrix_folder = "optimal_matrices/consumption/Nr25_Nc25/"
G_matrix_save_path = "optimal_matrices/consumption/Nr25_Nc25/RandTrix_Nr25_Nc25_Nest0.4_Conn0.1232.txt"
A_matrix_folder = "data_output/"
A_matrix_suffix = "_alpha0=1_intra_specific_syntrophy=not_allowed_verbose-level=1_gamma0=1_energy"
Npoints = 2000


columns=["NR", "NS", "connG", "nestG", "E", "connA", "nestA"]
df = pd.DataFrame(columns=columns)

G_matrix_name = G_matrix_save_path[len(G_matrix_folder):]
A_matrix_name = G_matrix_name[:-4]+"_optimal_alpha.txt"
[NR, NS, nestG, connG] = cf.extract_numbers_from_string(G_matrix_name)
nestG = cf.closest_element_in_list(nestG, cf.all_nestedness)
connG = cf.closest_element_in_list(connG, cf.all_connectance)


A_data=np.loadtxt(A_matrix_folder+A_matrix_name+A_matrix_suffix)
energy = A_data[:, 0]
nestA = A_data[:, 1]
connA = A_data[:, 2]

df = df.append(pd.DataFrame([[NR, NS, connG, nestG, cf.asymptote(energy, Npoints), cf.asymptote(connA, Npoints), cf.asymptote(nestA, Npoints)]], columns=columns))
df.to_csv("data_output/optimized_matrices_characteristics.csv", index=None)
