import numpy as np
import matplotlib.pyplot as plt
import consumer_resource_data_analysis as cf
from consumer_resource_data_analysis import alpha_mode, label, alpha0, all_nestedness, all_connectance, alpha_mode_colours, nestedness_label, connectance_label
import os
import matplotlib.tri as tr
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit

import copy

filename = 'local_dynamical_stability/local_dynamical_stability_NR25_NS25_100_points_full_rank_opt_consumption_mat_NR25_NS25'
optimal_LRI_folder='optimal_LRI'
consumption_matrix_folder='optimal_matrices/consumption/Nr25_Nc25'

cmap = plt.cm.get_cmap('jet_r')
colors = [cmap(i/10) for i in range(len(alpha0))]

# FIRST LOAD DATA
# local_dynamical_stability region[alpha_mode][alpha0][connectance][nestedness][gamma0][S0] contains the local_dynamical_stability of said point
# cf.filter_data(alpha_mode, alpha0, filename, optimal_LRI_folder, consumption_matrix_folder)
local_dynamical_stability_region=cf.load_data_region(alpha_mode,alpha0, filename,optimal_LRI_folder)
alpha0=np.array(alpha0)

# then we plot the heat map of the dynamical stability probability
index_matrix = 25
index_alpha_mode = 0
index_alpha0 = [0,1,9]


pmin = 0.8
levels = np.linspace(pmin, 1., 100)

fig, axs = plt.subplots(1, 3, sharey=True, sharex=True, figsize=(10.5,4.5))

for i in range(len(index_alpha0)):
    data = local_dynamical_stability_region[index_alpha_mode, index_alpha0[i], index_matrix]
    NR = int(data[0])
    NS = int(data[1])
    nest = data[2]
    conn = data[3]
    gamma0 = data[4::3]
    S0 = data[5::3]
    proba_lds = data[6::3]

    Npoints = int(np.sqrt(len(proba_lds)))
    proba_lds = np.transpose(np.reshape(proba_lds,(Npoints,Npoints)))
    gamma0 = np.linspace(np.min(gamma0), np.max(gamma0), Npoints)
    S0 = np.linspace(np.min(S0), np.max(S0), Npoints)
    ax = axs[i]
    ax.set_aspect('equal')
    im = ax.contourf(gamma0, S0, proba_lds, levels=levels, cmap=cmap)

    ax.set_xlabel(r'$\gamma_0$')
    ax.set_xlim(0.01, 1)
    ax.set_ylim(0.01, 1)
    ax.set_title(r'$\alpha_0='+str(alpha0[index_alpha0[i]])+'$')

    axs[i].set_xticks([0.01, 0.5, 1])
    axs[i].set_xticklabels([0.01, 0.5, 1])

axs[0].set_ylabel(r'$S_0$')
axs[0].set_yticks([0.01,0.5, 1])
axs[0].set_yticklabels([0.01,0.5, 1])
fig.subplots_adjust(bottom=0.2, top=0.95)

cbar_ax = fig.add_axes([0.125, 0.15, 0.75, 0.02])
cbar=fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
cbar.set_ticks(np.linspace(pmin, 1, 5))
cbar.set_ticklabels(np.round(np.linspace(pmin, 1, 5), decimals=3))
cbar.set_label(r'$\mathcal{D}_L\left(\gamma_0, S_0, G\right)$')
title = r'$N_R='+str(NR)+', N_S='+str(NS)+', \eta_G='+str(round(nest,2))+', \kappa_G='+str(round(conn,2))+'$'
fig.suptitle(title)
savename ='NR'+str(NR)+'_NS'+str(NS)+'_Nest'+str(nest)+'_Conn'+str(conn)
fig.savefig('plots/probability_dynamical_stability_'+savename+'.pdf')
