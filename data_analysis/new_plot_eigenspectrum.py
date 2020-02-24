import numpy as np
import matplotlib.pyplot as plt
from common_functions import remove_strings_from_file

file_folder='data_output'
file_name='test'
save_folder='plots'
save_name='eigenvalues_Nest0.15_Conn01808_sigma_spread'

file_to_load = file_folder+'/'+file_name
remove_strings_from_file('matrices/Nr25_Nc25', file_to_load)
data = np.loadtxt(file_to_load+'_filtered.out', dtype=complex)

NR=data[0]
NS=data[1]
nestedness=data[2]
connectance=data[3]
eigenvalues=np.sort_complex(data[4:])


size= []
color = []
for j in range(len(eigenvalues)):
    if np.real(eigenvalues[j]) > 1e-16:
        color.append('red')
        size.append(10)
    else:
        color.append('blue')
        size.append(0.2)



fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
ax1.scatter(np.real(eigenvalues), np.imag(eigenvalues), s=size, color=color)

fig1.tight_layout()
fig1.savefig(save_folder+'/'+save_name+".pdf")
