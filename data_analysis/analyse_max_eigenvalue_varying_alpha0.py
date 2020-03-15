import numpy as np
import matplotlib.pyplot as plt
import sys

file_to_load = sys.argv[1]
output = sys.argv[2]
matrices_folder = 'optimal_matrices/consumption/Nr25_Nc25'
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

marginally_stable=[]
always_unstable=[]
always_stable=[]
transition_observed=[]

# we write the maximum stability and feasability for each matrix
file=open(output, 'w')
# for each matrix, we estimate with a linear fit when max(lambda)(alpha0)=0
for i in range(len(data)):
    file.write(str(int(NR[i]))+' '+str(int(NS[i]))+' '+str(nestedness[i])+' '+str(connectance[i])+' ')
    string = "Matrix with nestedness "+str(nestedness[i])+ " and connectance " + str(connectance[i])
    if max_l[i,0]*max_l[i,-1]>0:
        if max_l[i,-1] < -ZERO:
            always_stable.append(i)
            string+=" is always dynamically stable"
            file.write('stable')
        elif max_l[i,-1] > ZERO:
            always_unstable.append(i)
            string+=" is always dynamically unstable"
            file.write('unstable')
        else:
            marginally_stable.append(i)
            string+=" is marginally stable"
    else:
        transition_observed.append(i)
        # gives the index of the last negative element
        index=0
        while(max_l[i][index]<0):
            index+=1
        max_dyn_stable_a0=0.5*(alpha0[i][index]+alpha0[i][index+1])
        file.write(str(max_dyn_stable_a0))
        string += ", maximally dynamically stable a0 is " + str(max_dyn_stable_a0)
    file.write(' '+str(alpha0[i,-1]))
    file.write('\n')
file.close()
