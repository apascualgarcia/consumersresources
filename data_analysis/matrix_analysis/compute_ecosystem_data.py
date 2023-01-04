import numpy as np
from consumer_resource_data_analysis import mat_connectance, mat_nestedness, eco_energy, all_nestedness, all_connectance, closest_element_in_list
import pandas as pd


def add_eco_data_to_df(df, G_path, A_path, NR, NS, columns):
    G = np.loadtxt(G_path)
    A = np.loadtxt(A_path)
    to_append = pd.DataFrame([[NR, NS, closest_element_in_list(mat_connectance(G), all_connectance), closest_element_in_list(mat_nestedness(G), all_nestedness), eco_energy(A,G), mat_connectance(A), mat_nestedness(A)]],
                columns = columns)
    df = df.append(to_append)

    return df



G_matrix_folder = "optimal_matrices/consumption/Nr25_Nc25/"
G_matrix_list = "matrix_list/full_rank_opt_consumption_mat_NR25_NS25.in"
A_matrix_folder = "optimal_matrices/syntrophy/Nr25_Nc25/9Jul21//"
A_matrix_suffix = "_optimize_matrices_alpha0=1_intra_specific_syntrophy=allowed_gamm_no_iss.txt"
columns=["NR", "NS", "connG", "nestG", "E", "connA", "nestA"]
save_file = "data_output/optimized_matrix_data_alpha0=1_intra_specific_syntrophy=allowed_no_iss.csv"
df = pd.DataFrame(columns=columns)

NR=25
NS=25

for G_matrix_save_path in np.loadtxt(G_matrix_list, dtype='U'):
    G_matrix_name = G_matrix_save_path[len(G_matrix_folder):]
    A_matrix_name = G_matrix_name[:-4]+"_optimal_alpha.txt"
    A_matrix_path = A_matrix_folder+A_matrix_name+A_matrix_suffix

    df = add_eco_data_to_df(df, G_matrix_save_path, A_matrix_path, NR, NS, columns)

df.to_csv(save_file, index=None, sep="\t")
