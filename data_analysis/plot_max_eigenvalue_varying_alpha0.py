import numpy as np
import matplotlib.pyplot as plt
import sys

file_name = 'max_eigenvalues_variable_syntrophy_full_rank_opt_consumption_mat_NR25_NS25_optimal_matrix_optimal_LRI'
matrices_folder = 'optimal_matrices/consumption/Nr25_Nc25'
ZERO = 1e-15

index_mat = 1


data_folder='data_output'


file_to_load = data_folder+'/'+file_name
def remove_strings_from_file(filename):
    file = open(filename + '.out', "r")
    metadata = []
    for x in file:
        name = x.replace(matrices_folder + '/RandTrix_Nr', '')
        name = name.replace('_Nc', ' ')
        name = name.replace('_Nest', ' ')
        name = name.replace('_Conn', ' ')
        name = name.replace('.txt', '')
        metadata.append(name.rstrip())
    file.close()
    f = open(filename + '_filtered.out', 'w')
    for a in metadata:  # python will convert \n to os.linesep
        f.write(a + '\n')
    f.close()

remove_strings_from_file(file_to_load)
data = np.loadtxt(file_to_load+'_filtered.out')

NR=data[:,0]
NS=data[:,1]
nestedness=data[:,2]
connectance=data[:,3]
alpha0=data[:, 4::2]
max_l = data[:, 5::2]


color=[]
for i in range(len(alpha0[index_mat])):
    if max_l[index_mat][i] > 0:
        color.append('red')
    else:
        color.append('blue')


fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
ax1.plot(alpha0[index_mat], max_l[index_mat])
ax1.set_yscale('symlog')

ax1.set_title('')

plt.show()
