import numpy as np
import matplotlib.pyplot as plt

file_folder='data_output'
file_name = 'largest_eigenvalues_varying_alpha0'
matrices_folder = 'matrices/Nr25_Nc25'

def remove_strings_from_file(filename):
    file = open(filename + '.out', "r")
    metadata = []
    for x in file:
        name = x.replace(matrices_folder + '/RandTrix_Nr', '')
        name = name.replace('_Nc', ' ')
        name = name.replace('_Nest', ' ')
        name = name.replace('_Conn', ' ')
        name = name.replace('.txt', '')
        metadata.append(name)
    file.close()
    f = open(filename + '_filtered.out', 'w')
    for a in metadata:  # python will convert \n to os.linesep
        f.write(a + '\n')
    f.close()

file_to_load = file_folder+'/'+file_name
remove_strings_from_file(file_to_load)

data = np.loadtxt(file_to_load+'_filtered.out')

NR=data[:,0]
NS=data[:,1]
nestedness=data[:,2]
connectance=data[:,3]
alpha0=data[:, 4::2]
max_l = data[:, 5::2]

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
ax1.plot(alpha0[2], max_l[2])

fig1.tight_layout()
plt.show()
