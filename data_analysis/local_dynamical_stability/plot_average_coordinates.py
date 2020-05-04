import numpy as np
import matplotlib.pyplot as plt
import consumer_resource_data_analysis as cf
from consumer_resource_data_analysis import alpha_mode, label, alpha0, all_nestedness, all_connectance, alpha_mode_colours, nestedness_label, connectance_label
import os
import matplotlib.tri as tr
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit

import copy

filename = 'local_dynamical_stability/leo_file_local_dynamical_stability_NR25_NS25_100_points_full_rank_opt_consumption_mat_NR25_NS25'
optimal_LRI_folder='optimal_LRI_Nr25_Nc25'
consumption_matrix_folder='optimal_matrices/consumption/Nr25_Nc25'
matrix_set='S_{25}'

cmap = plt.cm.get_cmap('jet_r')
colors = [cmap(i/10) for i in range(len(alpha0))]

# FIRST LOAD DATA
# local_dynamical_stability region[alpha_mode][alpha0][connectance][nestedness][gamma0][S0] contains the local_dynamical_stability of said point
cf.filter_data(alpha_mode, alpha0, filename, optimal_LRI_folder, consumption_matrix_folder)
local_dynamical_stability_region=cf.load_data_region(alpha_mode,alpha0, filename,optimal_LRI_folder)
alpha0=np.array(alpha0)

fig=plt.figure()
ax = fig.add_subplot(111)
for i in range(len(alpha_mode)):
    print('alpha mode : ', alpha_mode[i])
    av_gamma0=[]
    for k in range(len(alpha0)):
        print('alpha0 = ', alpha0[k])
        av_gamma0.append(cf.average_gamma0_all_matrices(local_dynamical_stability_region, i, k))
    ax.plot(alpha0, av_gamma0, color=alpha_mode_colours[i], label=label[i])
ax.set_xlabel(r'$\alpha_0$')
ax.set_ylabel(r'$\langle \gamma_0 \rangle_D (\alpha_0)$')
ax.legend()
fig.tight_layout()
fig.savefig('plots/local_dynamical_stability_average_gamma0.pdf')

fig=plt.figure()
ax = fig.add_subplot(111)
for i in range(len(alpha_mode)):
    av_S0=[cf.average_S0_all_matrices(local_dynamical_stability_region, i, k) for k in range(len(alpha0))]
    ax.plot(alpha0, av_S0, color=alpha_mode_colours[i], label=label[i])
ax.set_xlabel(r'$\alpha_0$')
ax.set_ylabel(r'$\langle S_0 \rangle_D (\alpha_0)$')
ax.legend()
fig.tight_layout()
fig.savefig('plots/local_dynamical_stability_average_S0.pdf')
