import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d


data_folder="data_output"
file="random_matrix_maximum_feasability"


filename = data_folder+"/"+file+".out"
data = np.loadtxt(filename)

fig1=plt.figure(1)
ax1 = fig1.add_subplot(111, projection='3d')
ax1.plot(data[:,0], data[:,1], data[:,2])
ax1.set_xlabel(r'$S_0$')
ax1.set_ylabel(r'$\gamma_0$')
ax1.set_zlabel(r'Max feasible $\alpha_0$')

plt.show()
