import numpy as np
import os
from consumer_resource_data_analysis import alpha_mode

strstab_folder='structural_stability'
start_string='strstab'
syntrophy_list=['common_max_syntrophy','own_max_syntrophy', 'no_syntrophy']
save_output='data_all_structural_stability'
header='organization of data : common max syntrophy (fully, no copr, optimal, random struct), own max syntrophy (F, NC, O, R) and no syntrophy'
files_list=[a for a in os.listdir(strstab_folder) if a.startswith(start_string)]
matrix_list=[a[len(start_string)+1:] for a in files_list]
files_w_folder=[strstab_folder+'/'+a for a in files_list]

for i in range(len(matrix_list)):
    a = matrix_list[i].split("_",2)
    b = a[0]+'_'+a[1]
    matrix_list[i] = b
matrix_name = set(matrix_list)
matrix_data=[]
for mat in matrix_name:
    # take all 9 files relevant to the matrix
    mat_files = [a for a in files_w_folder if mat in a]
    mat_name = np.loadtxt(mat_files[0], usecols=0, dtype='U')
    line = [str(mat_name)]
    for syn in syntrophy_list:
        # separate them by syntrophy type
        syn_files=[a for a in mat_files if syn in a]
        if syn!='no_syntrophy':
            for alpha in alpha_mode:
                relevant_file=[a for a in syn_files if alpha in a]
                to_add = np.loadtxt(relevant_file[0], usecols=(1,2))
                if len(to_add)==0:
                    to_add = np.array(['NaN', 'NaN'])
                line.append(to_add[0])
                line.append(to_add[1])

        else:
            relevant_file=syn_files
            to_add = np.loadtxt(relevant_file[0], usecols=(1,2))
            line.append(to_add[0])
            line.append(to_add[1])
    matrix_data.append(line)
np.savetxt(strstab_folder+'/'+save_output+'.out',matrix_data, fmt='%s', header=header)
