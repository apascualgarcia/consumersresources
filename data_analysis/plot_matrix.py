import matplotlib.pyplot as plt
import numpy as np

matrices_folder='matrices/Nr25_Nc25'
matrix_name = 'RandTrix_Nr25_Nc25_Nest0.2_Conn0.1312'
#matrix_name = 'optimal_alpha_'+matrix_name

mat = np.loadtxt(matrices_folder+'/'+matrix_name+'.txt')

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
ax1.matshow(mat)

plt.show()
