import numpy as np
import glob
import consumer_resource_data_analysis as cf
import pandas as pd

########### NON CUSTOMIZABLE PART ################
def append_one_matrix_data(df, mat_files):
    E = []
    K_A =[]
    K_G =[]
    ETA_A=[]
    ETA_G=[]

    for file in mat_files:
        data = np.loadtxt(file)
        if data.size > 0:
            energy = data[0]
            nestA = data[1]
            connA = data[2]
            nestG = data[3]
            connG = data[4]

            E.append(energy)
            K_A.append(connA)
            K_G.append(connG)
            ETA_A.append(nestA)
            ETA_G.append(nestG)
    data ={'NR':25,'NS':25,'<connG>':np.mean(K_G),'<connA>':np.mean(K_A),'<E>': np.mean(E),
        '<nestG>':np.mean(ETA_G),'<nestA>':np.mean(ETA_A),'connG std':np.std(K_G, ddof=1),
        'connA std':np.std(K_A, ddof=0),'E std':np.std(E, ddof=0),'nestG std':np.std(ETA_G, ddof=0),
        'nestA std':np.std(ETA_A, ddof=0)}
    print("For kg=", np.mean(K_G), "E=", E)
    df = df.append(data, sort=False, ignore_index=True)
    return df


mat_list = "matrix_list/full_rank_opt_consumption_mat_NR25_NS25.in"
mat_folder = "optimal_matrices/consumption/Nr25_Nc25"
data_folder = "optimal_matrices/consumption/Nr25_Nc25"
suffix = "_energy"
save_file = "data_output/optimized_matrix_data_all.csv"
columns=['NR', 'NS', '<connG>', '<nestG>', '<E>', '<connA>', '<nestA>', 'connG std', 'nestG std', 'E std', 'connA std', 'nestA std']


df = pd.DataFrame(data=None, columns=columns)
for G_mat_path in np.loadtxt(mat_list, dtype='U'):
    G_name = data_folder+G_mat_path[len(mat_folder):-4]+"*"+suffix
    files = glob.glob(G_name)
    df = append_one_matrix_data(df, files)
df.to_csv(save_file, index=None)
