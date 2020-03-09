import matplotlib.pyplot as plt
import numpy as np

gamma_folder_path='optimal_matrices/consumption/Nr25_Nc25/RandTrix_Nr25_Nc25_Nest0.1_Conn0.1296.txt'
alpha_folder_path='optimal_matrices/optimal_LRI/RandTrix_Nr25_Nc25_Nest0.1_Conn0.1296_optimal_alpha.txt'

alpha0=0.75
gamma0=0.1

G = np.loadtxt(gamma_folder_path)
A = np.loadtxt(alpha_folder_path)

O = A@G
C = np.transpose(G)@G

fig1=plt.figure(1)
ax1=fig1.add_subplot(111)
im=ax1.imshow(alpha0*O-gamma0*C)
cbar=fig1.colorbar(im)

plt.show()
