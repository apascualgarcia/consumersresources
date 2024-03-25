import numpy as np
from os import listdir
from os.path import isfile, join
import re
from datetime import datetime


origin_folder="./matrices/Nr25_Nc25/syntrophy/optimized_alpha01/Binary"
destination_folder="./matrices/Nr25_Nc25/syntrophy/optimized_alpha01/Metamatrices"

matrix_list = "./matrix_list/Nr25_Nc25/consumption/full_rank_opt_consumption_mat_NR25_NS25.in"
seeds = range(50)

matrices = np.genfromtxt(matrix_list, dtype='str')
for mat in matrices:
    start_name = mat.find("RandTrix")
    mat_name = mat[start_name:-4:]

    act_seeds = []
    opt_mats = []

    for seed in seeds:
        root_folder = origin_folder+"/seed_"+str(seed)
        opt_matrices = listdir(root_folder)
        opt_mat = [root_folder+"/"+str(x) for x in opt_matrices if (mat_name in x and x[-6::]!="energy")]

        if len(opt_mat)==1:
            act_seeds.append(seed)
            opt_mats.append(np.loadtxt(opt_mat[0]))

    if len(opt_mats) > 0:
        #### on moyenne pour obtenir une métamatrice ####
        meta_mat = opt_mats[0]
        for i in range(len(opt_mats)-1):
            meta_mat+=opt_mats[i]
        meta_mat = np.divide(meta_mat, len(opt_mats))
        
        now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        header="created "+now+" with binary optimized matrices with seeds "+str(act_seeds)+" for matrix "+mat_name
        np.savetxt(destination_folder+"/"+mat_name+"_optimal_alpha.txt",meta_mat,fmt="%.5f", header=header)
    else:
        print(mat_name)
        