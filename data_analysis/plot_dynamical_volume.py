import numpy as np
import matplotlib.pyplot as plt
from common_functions import remove_strings_from_file
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import sys
import copy
import matplotlib.tri as tr
from scipy.optimize import curve_fit
from matplotlib import rcParams


np.set_printoptions(threshold=sys.maxsize)

data_folder='./data_output'
plot_folder='./plots'
matrices_folder='optimal_matrices/consumption/Nr25_Nc25'
filename='local_stability_NR25_NS25'
mat_index=4
mat_index_2=23
exp_coeff=[0.2,0.25,0.111111,0.111111,0.125,0.0833333,0.0833333,0.0833333,0.1,0.111111,0.0769231,0.0666667,0.0666667,0.0769231,0.0769231,0.0714286,0.0769231,0.0588235,0.0588235,0.0625,0.0666667,0.0625,0.0666667,0.0555556,0.0555556,0.0526316]

filename = data_folder+'/'+filename
remove_strings_from_file(matrices_folder, filename)
data=np.loadtxt(filename+'_filtered.out')

def func_fit(x, a,b):
    return a/x+b


NR=data[:,0]
NS=data[:,1]
nestedness=data[:,2]
connectance=data[:,3]
gamma0=data[:, 4::3]
S0=data[:, 5::3]
stability=data[:, 6::3]

# plot typical matrix stability volume
title_axis=r'$N_R='+str(int(NR[mat_index]))
title_axis+=r', N_S='+str(int(NS[mat_index]))
title_axis+=r', \eta='+str(nestedness[mat_index])
title_axis+=r', \kappa='+str(connectance[mat_index])
title_axis+='$'

title_plot=r'$N_R='+str(int(NR[mat_index]))
title_plot+=r', N_S='+str(int(NS[mat_index]))+'$'

new_y_coord = np.array([gamma0[mat_index][i]*(S0[mat_index][i]+0.0046)/0.043 for i in range(len(gamma0[mat_index]))])

cmap = cm.get_cmap('bwr_r')
fig1 =plt.figure(1)
ax1 = fig1.add_subplot(111, aspect='equal')
#im = ax1.tricontourf(gamma0[mat_index], S0[mat_index], stability[mat_index], cmap=cmap, levels=np.linspace(0., 1.001, 10))
im=ax1.scatter(gamma0[mat_index], new_y_coord, c=stability[mat_index])

ax1.set_xlim(0., 1.)
ax1.set_ylim(0., 1.)
ax1.set_xlabel(r'$\gamma_0$')
ax1.set_ylabel(r'$\gamma_0S_0$')
title_axis=r'$\eta='+str(nestedness[mat_index])
title_axis+=r', \kappa='+str(connectance[mat_index])+'$'
ax1.set_title(title_axis)


print(stability)

#fig1.subplots_adjust(right=3)
cbar=fig1.colorbar(im, ax=[ax1], orientation='horizontal', shrink=0.8, pad=-0.4)
cbar.set_label(r'$\mathcal{F}(\gamma_0, S_0, G)$')
fig1.suptitle(title_plot)
fig1.tight_layout()

plt.show()
