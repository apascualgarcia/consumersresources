import numpy as np
import matplotlib.pyplot as plt
import sys

file_to_load = sys.argv[1]
matrices_folder = './optimal_matrices/Nr25_Nc25'
ZERO = 1e-15

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

remove_strings_from_file(file_to_load)
data = np.loadtxt(file_to_load+'_filtered.out')

NR=data[:,0]
NS=data[:,1]
nestedness=data[:,2]
connectance=data[:,3]
alpha0=data[:, 4::2]
max_l = data[:, 5::2]

marginally_stable=[]
always_unstable=[]
always_stable=[]
transition_observed=[]

# for each matrix, we estimate with a linear fit when max(lambda)(alpha0)=0
for i in range(len(data)):
    if max_l[i,0]*max_l[i,-1]>0:
        if max_l[i,-1] < -ZERO:
            always_stable.append(i)
        elif max_l[i,-1] > ZERO:
            always_unstable.append(i)
        else:
            marginally_stable.append(i)
    else:
        transition_observed.append(i)
        linear_fit=np.polyfit(alpha0, max_l, 1)
        max_dyn_stable_a0=-linear_fit[1]/linear_fit[0]
        print("For this matrix, maximal dynamically stable a0 is ", max_dyn_stable_a0)
