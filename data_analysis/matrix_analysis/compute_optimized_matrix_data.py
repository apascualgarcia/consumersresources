import numpy as np
import pandas as pd
import consumer_resource_data_analysis as cf

 ##### NON CUSTOMIZABLE PART ####################
def append_mat_data_to_df(df, G_matrix_name, G_matrix_folder, A_matrix_path, columns, Npoints):

    [NR, NS, nestG, connG] = cf.extract_numbers_from_string(G_matrix_name)
    nestG = cf.closest_element_in_list(nestG, cf.all_nestedness)
    connG = cf.closest_element_in_list(connG, cf.all_connectance)

    A_data=np.loadtxt(A_matrix_path)
    energy = A_data[:, 0]
    nestA = A_data[:, 1]
    connA = A_data[:, 2]

    df = df.append(pd.DataFrame([[int(NR), int(NS), connG, nestG, cf.asymptote(energy, Npoints), cf.asymptote(connA, Npoints), cf.asymptote(nestA, Npoints)]], columns=columns))
    return df

def append_mat_from_Gnname(df, G_matrix_folder, G_matrix_save_path, A_matrix_folder, A_matrix_suffix, columns, Npoints):
    G_matrix_name = G_matrix_save_path[len(G_matrix_folder):]
    A_matrix_name = G_matrix_name[:-4]+"_optimal_alpha.txt"
    A_matrix_path = A_matrix_folder+A_matrix_name+A_matrix_suffix

    df = append_mat_data_to_df(df, G_matrix_name, G_matrix_folder, A_matrix_path, columns, Npoints)
    return df
###### END OF NON CUSTOMIZABLE PART ##########


G_matrix_folder = "optimal_matrices/consumption/Nr25_Nc25/"
G_matrix_list = "matrix_list/full_rank_opt_consumption_mat_NR25_NS25.in"
A_matrix_folder = "data_output/"
A_matrix_suffix = "_optimize_matrices_alpha0=1_intra_specific_syntrophy=allowed_gamma0=1_energy_5millionssteps"
save_file = "data_output/optimized_matrix_data_alpha0=1_intra_specific_syntrophy=allowed_5millionssteps.csv"
Npoints = 2000
columns=["NR", "NS", "connG", "nestG", "E", "connA", "nestA"]
df = pd.DataFrame(columns=columns)


for G_matrix_save_path in np.loadtxt(G_matrix_list, dtype='U'):
    try:
        df = append_mat_from_Gnname(df, G_matrix_folder, G_matrix_save_path, A_matrix_folder, A_matrix_suffix, columns, Npoints)
    except:
        continue
df.to_csv(save_file, index=None, sep="\t")
