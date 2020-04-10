import numpy as np
import os

origin_folder='optimal_matrices/consumption/Nr25_Nc25'
output_folder='optimal_matrices/consumption/list_representation'

files = [f for f in os.listdir(origin_folder) if not f.startswith('.')]
for f in files:
    mat = np.loadtxt(origin_folder+'/'+f)
    list_rep = []
    for i in range(len(mat)):
        for mu in range(len(mat[0])):
            if(mat[i,mu]==1):
                string='Species'+str(i+1)+' Resource'+str(mu+1)
                list_rep.append(string)
    list_rep = np.transpose(np.array(list_rep))
    np.savetxt(output_folder+'/'+f, list_rep, fmt='%s')
