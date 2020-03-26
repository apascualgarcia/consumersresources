import numpy as np
import matplotlib.pyplot as plt
import common_functions as cf
import os
import matplotlib.tri as tr
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit
from matplotlib.colors import LogNorm

import copy

alpha_mode=['random_structure', 'no_release_when_eat', 'optimal_matrix']
label=['fully connected', 'no release when eat', 'optimal LRI']
filename = 'largest_eigenvalue/largest_eigenvalue_NR25_NS25_full_rank_opt_consumption_mat_NR25_NS25'
alpha0=[0, 1.3e-3, 2.6e-3, 3.9e-3, 5.2e-3, 6.5e-3, 7.8e-3, 9.1e-3, 1.04e-2, 1.4e-2]
all_nestedness=[0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]
all_connectance=[0.08, 0.13, 0.18, 0.23, 0.28, 0.33, 0.38, 0.43]
Npoints = 900


cmap = plt.cm.get_cmap('jet_r')
colors = [cmap(i/10) for i in range(len(alpha0))]

# largest_eigenvalue region[alpha_mode][alpha0][connectance][nestedness][gamma0][S0] contains the largest_eigenvalue of said point
largest_eigenvalue_region = []
for al_mo in alpha_mode:
    local_vector=[]
    for a in alpha0:
        # file = filename+'_'+al_mo+'_optimal_LRI_alpha0='+str(a)
        # cf.remove_strings_from_file('optimal_matrices/consumption/Nr25_Nc25', file)
        file = filename+'_'+al_mo+'_optimal_LRI_alpha0='+str(a)+'_filtered.out'
        local_data=np.loadtxt(file, dtype=complex)
        local_vector.append(local_data)
    largest_eigenvalue_region.append(local_vector)

largest_eigenvalue_region = np.array(largest_eigenvalue_region)
connectance = np.real(largest_eigenvalue_region[0,0][:,3])
nestedness = np.real(largest_eigenvalue_region[0,0][:,2])

alpha0=np.array(alpha0)

# plot largest eigenvalue plot at no syntrophy for all matrices
# for k in range(len(largest_eigenvalue_region[0,0])):
#     for j in range(len(largest_eigenvalue_region[0])):
#         data = np.ma.masked_invalid(largest_eigenvalue_region[:,j,k])
#         if(not(data[:,6::3].mask.all())):
#             fig, axs = plt.subplots(1, 3, sharey=True, sharex=True, figsize=(10.5,4.5))
#             max_ev=-1000
#             min_ev=1000
#             # first find the largest and smallest eigenvalues for the three regimes
#             for i in range(len(largest_eigenvalue_region)):
#                 largest_ev=data[i,6::3]
#                 if(np.min(np.real(largest_ev)) < min_ev):
#                     min_ev = np.min(np.real(largest_ev))
#                 if(np.max(np.real(largest_ev)) > max_ev):
#                     max_ev = np.max(np.real(largest_ev))
#
#             for i in range(len(largest_eigenvalue_region)):
#                 gamma0=np.real(data[i,4::3])
#                 S0=np.real(data[i,5::3])
#                 NR=np.real(data[i,0])
#                 NS=np.real(data[i,1])
#                 nestedness=np.real(data[i,2])
#                 connectance=np.real(data[i,3])
#                 largest_ev=data[i,6::3]
#
#
#                 axs[i].set_aspect('equal')
#                 im=axs[i].scatter(gamma0, S0, c=np.real(largest_ev), s=25, marker='s', vmin=min_ev, vmax=max_ev, cmap='jet')
#                 axs[i].set_xlabel(r'$\gamma_0$')
#                 axs[i].set_xticks([0, 0.5, 1])
#                 axs[i].set_xticklabels([0, 0.5, 1])
#                 axs[i].set_title(label[i])
#
#                 axs[i].set_xlim(0.,1.)
#                 axs[i].set_ylim(0.,1.)
#             axs[0].set_ylabel(r'$S_0$')
#             axs[0].set_yticks([0, 0.5, 1])
#             axs[0].set_yticklabels([0, 0.5, 1])
#             fig.subplots_adjust(bottom=0.2, top=0.95)
#             cbar_ax = fig.add_axes([0.125, 0.15, 0.75, 0.02])
#             cbar=fig.colorbar(im, cax=cbar_ax, orientation='horizontal', format='%.1e')
#             #cbar.set_label(r'$\alpha_0$')
#             title = r'Re$(\lambda_1)$ for $N_R='+str(int(NR))+', N_S='+str(int(NS))+', \kappa='+str(round(connectance,2))\
#                         +', \eta='+str(nestedness)+'$'
#             title +=r' at $\alpha_0='+str(alpha0[j])+'$'
#             fig.suptitle(title)
#             save_name='NR'+str(int(NR))+'_NS'+str(int(NS))+'_Nest'+str(nestedness)+'_Conn'+str(connectance)+'_alpha0='+str(alpha0[j])
#             #fig.savefig('plots/largest_eigenvalue_wt_wc_'+save_name+'.pdf')
#             plt.show()
#             plt.close()

critical_alpha0 = []
for j in range(len(largest_eigenvalue_region[0,0])):
    decline_volume = []

    # find largest alpha0 that has non-zero volume
    max_vol = 0
    i = 0
    continue_loop=True
    while continue_loop:
        data = np.ma.masked_invalid(largest_eigenvalue_region[0,i,j])
        if not(data[6::3].mask.all()):
            max_vol +=1
        else:
            continue_loop=False
        continue_loop= (continue_loop and (i < len(largest_eigenvalue_region[0])-1))
        i+=1

    volumes=[]
    for i in range(max_vol):
        data = np.ma.masked_invalid(largest_eigenvalue_region[0,i,j])
        volumes.append(data[6::3].count()/Npoints)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(alpha0[:max_vol], volumes)

    fit_function = cf.linear_function

    fit_volumes, popt, perror = cf.fit_data(fit_function, alpha0[:max_vol], volumes[:max_vol])
    alpha0_crit, alpha0_error = cf.zero_from_fit(fit_function, popt, perror)
    critical_alpha0.append(alpha0_crit)
    ax.plot(alpha0[:max_vol], fit_volumes, marker='None', linestyle='solid')


    ax.set_xlabel(r'$\alpha_0$')
    ax.set_ylabel(r'Vol($\mathcal{D}_{L,1}$) (normalized)')
    ax.set_title(r'$\alpha_0^*='+str(alpha0_crit)+'$')#'\pm '+str(alpha0_error)+'$')
    fig.tight_layout()
    #plt.show()
    plt.close()

# plot critical alpha0
fig = plt.figure()
ax = fig.add_subplot(111)
for nest in all_nestedness:
    indices = [int(i) for i in range(len(nestedness)) if cf.closest_element_in_list(nestedness[i], all_nestedness)==nest]
    sorted_indices=np.array([indices[a] for a in np.argsort(connectance[indices])])
    connectance_to_plot=[connectance[i] for i in sorted_indices]
    critical_alpha0_to_plot=[critical_alpha0[i] for i in sorted_indices]
    ax.plot(connectance_to_plot, critical_alpha0_to_plot, label=r'$\eta\approx'+str(nest)+'$')
ax.set_xlabel(r'Connectance $\kappa$')
ax.set_ylabel(r'$\alpha_0^*$')
# ax.set_title(label)
ax.legend(bbox_to_anchor=(1.0, 1.0))
fig.tight_layout()
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)
for conn in all_connectance:
    indices = [int(i) for i in range(len(connectance)) if cf.closest_element_in_list(connectance[i], all_connectance)==conn]
    sorted_indices=np.array([indices[a] for a in np.argsort(nestedness[indices])])
    nestedness_to_plot=[nestedness[i] for i in sorted_indices]
    critical_alpha0_to_plot=[critical_alpha0[i] for i in sorted_indices]
    ax.plot(nestedness_to_plot, critical_alpha0_to_plot, label=r'$\kappa\approx'+str(conn)+'$')
ax.set_xlabel(r'Ecological overlap $\eta$')
ax.set_ylabel(r'$\alpha_0^*$')
# ax.set_title(label[k])
ax.legend(bbox_to_anchor=(1.0, 1.0))

fig.tight_layout()
plt.show()
