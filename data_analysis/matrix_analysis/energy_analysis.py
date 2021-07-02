import common_features.mpl_params
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

matrix_data_path = "test_matrices/RandTrix_Nr25_Nc25_Nest0.55_Conn0.1248_optimal_alpha.txt_energy"

data = np.loadtxt(matrix_data_path)
energy = data[:,0]
nest = data[:,1]
conn = data[:, 2]
T = data[:, 3]

fig1 = plt.figure(1)
plt.plot(energy, 'o', markersize=0.5)
plt.plot(moving_average(energy, 15000), marker='', linestyle='solid', linewidth=2)
plt.ylabel(r'$E$(steps)')
plt.xlabel(r'steps')
plt.tight_layout()


fig2 = plt.figure(2)
plt.plot(nest, 'o', markersize=0.5)
plt.plot(moving_average(nest, 15000), marker='', linestyle='solid', linewidth=2)
plt.xlabel(r'steps')
plt.ylabel(r'$\eta_A$(steps)')
plt.tight_layout()

fig3 = plt.figure(3)
plt.plot(conn, 'o', markersize=0.5)
plt.plot(moving_average(conn, 15000), marker='', linestyle='solid', linewidth=2)
plt.xlabel(r'steps')
plt.ylabel(r'$\kappa_A$(steps)')
plt.tight_layout()

fig4 = plt.figure(4)
plt.plot(T, 'o', markersize=0.5)
plt.yscale('log')
plt.xlabel(r'steps')
plt.ylabel(r'$T$(steps)')


plt.tight_layout()
plt.show()
