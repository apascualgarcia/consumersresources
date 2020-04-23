import numpy as np
from consumer_resource_data_analysis import remove_strings_from_file
import os

matrices_folder='optimal_matrices/consumption/Nr50_Nc25'
filename='structural_stability/critical_dynamical_syntrophies_NR50_NS25_backup'

alpha_mode=['fully_connected', 'no_release_when_eat', 'optimal_matrix', 'random_structure']


files_max_own_folder='config/structural_stability/maximal_own_syntrophies'
files_common_max_folder='config/structural_stability/common_max_syntrophies'
files_no_syntrophy_folder='config/structural_stability/no_syntrophy'

# files_max_own_folder='config/structural_stability/test'
# files_common_max_folder='config/structural_stability/test'
# files_no_syntrophy_folder='config/structural_stability/test'

remove_strings_from_file(matrices_folder, filename)
data = np.loadtxt(filename+'_filtered.out')

NR = [int(a) for a in data[:,0]]
NS = [int(a) for a in data[:,1]]
nest = data[:,2]
conn = data[:,3]

critical_alpha = np.genfromtxt(filename+'.out', usecols=[1,2,3,4])
matrix_name = np.genfromtxt(filename+'.out', usecols=[0], dtype=np.dtype('U'))

# generate files for maximal own syntrophies
# for i in range(len(matrix_name)):
#     for j in range(len(alpha_mode)):
#         # first create appropriate file
#         filename='strstab_Nest'+str(nest[i])+'_Conn'+str(conn[i])+'_'+alpha_mode[j]
#         path_to_save='data_output/'+filename+'_own_max_syntrophy.out'
#         filename = files_max_own_folder+'/'+filename+'_own_max_syntrophy.in'
#         command='cp config/structural_stability/configuration_structural_stability.in '+filename
#         os.system(command)
#
#
#         # then add the extra lines needed
#         to_add = []
#         to_add.append('alpha0='+str(critical_alpha[i,j]))
#         to_add.append('alpha_mode='+alpha_mode[j])
#         to_add.append('path_to_food_matrix='+matrix_name[i])
#         to_add.append('path_to_save_file='+path_to_save)
#
#         for el in to_add :
#             command = 'echo "'+el+'" >> '+filename
#             os.system(command)

# generate files for commmon min syntrophy
common_alpha0 = np.min(critical_alpha)
print(common_alpha0)
# for i in range(len(matrix_name)):
#     for j in range(len(alpha_mode)):
#         # first create appropriate file
#         filename='strstab_Nest'+str(nest[i])+'_Conn'+str(conn[i])+'_'+alpha_mode[j]
#         path_to_save='data_output/'+filename+'_common_max_syntrophy.out'
#         filename = files_common_max_folder+'/'+filename+'_common_max_syntrophy.in'
#         command='cp config/structural_stability/configuration_structural_stability.in '+filename
#         os.system(command)
#
#
#         # then add the extra lines needed
#         to_add = []
#         to_add.append('alpha0='+str(common_alpha0))
#         to_add.append('alpha_mode='+alpha_mode[j])
#         to_add.append('path_to_food_matrix='+matrix_name[i])
#         to_add.append('path_to_save_file='+path_to_save)
#
#         for el in to_add :
#             command = 'echo "'+el+'" >> '+filename
#             os.system(command)
#
# # generate files for no syntrophy
# for i in range(len(matrix_name)):
#     # first create appropriate file
#     filename='strstab_Nest'+str(nest[i])+'_Conn'+str(conn[i])
#     path_to_save='data_output/'+filename+'_no_syntrophy.out'
#     filename = files_no_syntrophy_folder+'/'+filename+'_no_syntrophy.in'
#     command='cp config/structural_stability/configuration_structural_stability.in '+filename
#     os.system(command)
#
#
#     # then add the extra lines needed
#     to_add = []
#     to_add.append('alpha0='+str(0))
#     to_add.append('alpha_mode='+alpha_mode[0])
#     to_add.append('path_to_food_matrix='+matrix_name[i])
#     to_add.append('path_to_save_file='+path_to_save)
#
#     for el in to_add :
#         command = 'echo "'+el+'" >> '+filename
#         os.system(command)
