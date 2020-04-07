import numpy as np
import matplotlib.pyplot as plt
import consumer_resource_data_analysis as cf
from consumer_resource_data_analysis import alpha_mode, label, alpha0, all_nestedness, all_connectance
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
cf.filter_data(alpha_mode, alpha0, filename, optimal_LRI_folder, consumption_matrix_folder)
local_dynamical_stability_region=cf.load_data_region(alpha_mode,alpha0, filename,optimal_LRI_folder)
alpha0=np.array(alpha0)


# Plot common local_dynamical_stability volume
NR=int(local_dynamical_stability_region[0,0,0,0])
NS=int(local_dynamical_stability_region[0,0,0,1])
fig, axs, im, levels=cf.plot_common_region(local_dynamical_stability_region, alpha_mode, colors, label)
cbar=cf.add_colorbar_to_plot_levels(fig, im, levels, alpha0)
cbar.set_label(r'$\alpha_0$')
fig.suptitle(r'$\mathcal{D}_1^{S_M}$ for $N_R='+str(NR)+'$, $N_S='+str(NS)+'$')
plt.show()
#
# print("Passing to plots of local dynamical stability for each matrix ")
# # # now we produce a similar plot for all the matrices listed here
# for j in range(len(local_dynamical_stability_region[0,0])):
#     data = local_dynamical_stability_region[:,:,j]
#     NR=data[0,0,0]
#     NS=data[0,0,1]
#     nestedness=data[0,0,2]
#     connectance=data[0,0,3]
#
#     fig, axs, im, levels =cf.plot_levels(data, colors, label)
#     save_name='NR'+str(int(NR))+'_NS'+str(int(NS))+'_Nest'+str(nestedness)+'_Conn'+str(connectance)
#     # first save figure without color bar or title
#     fig.savefig('plots/local_dynamical_stability_wt_wc_region_'+save_name+'.pdf')
#
#     # add color bar to plot
#     cbar = cf.add_colorbar_to_plot_levels(fig, im, levels, alpha0)
#     cbar.set_label(r'$\alpha_0$')
#     fig.savefig('plots/local_dynamical_stability_wt_region_'+save_name+'.pdf')
#
#     # add suptitle to plot
#     fig.suptitle(r'Fully dynamically stable region $\mathcal{D}^G_{L,1}$ for $N_R='+str(int(NR))+', N_S='+str(int(NS))+', \kappa='+str(round(connectance,2))+', \eta='+str(nestedness)+'$')
#     fig.savefig('plots/local_dynamical_stability_region_'+save_name+'.pdf')
#
#     plt.close()
#
# decline =[]
# for j in range(len(local_dynamical_stability_region[0][0])):
#     fig, ax =  plt.subplots(1, 1)
#     local_decline=[]
#     for i in range(len(local_dynamical_stability_region)):
#         data=local_dynamical_stability_region[i,:,j]
#         NR=data[0,0]
#         NS=data[0,1]
#         nestedness=data[0,2]
#         connectance=data[0,3]
#         feas_volume = cf.shrink_volume_for_one_matrix(data)
#         ax.plot(alpha0, feas_volume, label=label[i],markersize=10, linewidth=2.5, markeredgewidth=3)
#         popt, pcov=curve_fit(cf.exp_fit, alpha0, feas_volume, maxfev=90000000)
#         # fitted_vol=[cf.exp_fit(x,*popt) for x in alpha0]
#         local_decline.append(popt[1])
#         #print(popt)
#         #ax.plot(alpha0, fitted_vol, linestyle='-', marker='None')
#     ax.set_xlabel(r'$\alpha_0$')
#     ax.set_ylabel(r'Vol$(\mathcal{D}^G_{L,1}(\alpha_0))$ (normalized)')
#     ax.set_ylim(-0.02,0.55)
#     ax.legend()
#     title=r'$N_R='+str(int(NR))+', N_S='+str(int(NS))+', \kappa='+str(round(connectance,2))+', \eta='+str(round(nestedness,2))+'$'
#     ax.set_title(title)
#     save_name='NR'+str(int(NR))+'_NS'+str(int(NS))+'_Nest'+str(nestedness)+'_Conn'+str(connectance)
#     fig.tight_layout()
#     print("Saving fig to "+'plots/size_local_dynamical_stability_region_'+save_name+'.pdf')
#     fig.savefig('plots/size_local_dynamical_stability_region_'+save_name+'.pdf')
#     decline.append(local_decline)
#     plt.close()
#     #plt.show()
# decline=np.array(decline)
# #
# # plot a slope curve vs nestedness and connectance
# decline = np.transpose(decline)
#
# connectance = local_dynamical_stability_region[0,0][:,3]
# nestedness = local_dynamical_stability_region[0,0][:,2]
#
# print("Passing to plots of exponential fit")
#
# for k in range(len(local_dynamical_stability_region)):
#     fig = plt.figure(k)
#     ax = fig.add_subplot(111)
#     for nest in all_nestedness:
#         indices = [i for i in range(len(nestedness)) if cf.closest_element_in_list(nestedness[i], all_nestedness)==nest]
#         ax.plot(connectance[indices],decline[k][indices], label=r'$\eta\approx'+str(nest)+'$')
#     ax.set_xlabel(r'Connectance $\kappa$')
#     ax.set_ylabel(r'Exponential coefficient')
#     ax.set_title(label[k])
#     ax.legend(bbox_to_anchor=(1.0, 1.0))
#     fig.tight_layout()
#     print("Saving fig to "+'plots/local_dynamical_stability_NR'+str(int(NR))+'_NS'+str(int(NS))+'_exp_fit_fixed_nestedness_'+label[k]+'.pdf')
#     fig.savefig('plots/local_dynamical_stability_NR'+str(int(NR))+'_NS'+str(int(NS))+'_exp_fit_fixed_nestedness_'+label[k]+'.pdf')
#     plt.close(k)
#
#     fig = plt.figure(k)
#     ax = fig.add_subplot(111)
#     for conn in all_connectance:
#         indices = [i for i in range(len(connectance)) if cf.closest_element_in_list(connectance[i], all_connectance)==conn]
#         sorted_indices=[ indices[a] for a in np.argsort(nestedness[indices])]
#         ax.plot(nestedness[sorted_indices],decline[k][sorted_indices], label=r'$\kappa\approx'+str(conn)+'$')
#     ax.set_xlabel(r'Nestedness $\eta$')
#     ax.set_ylabel(r'Exponential coefficient')
#     ax.set_title(label[k])
#     ax.legend(bbox_to_anchor=(1.0, 1.0))
#
#     fig.tight_layout()
#     print("Saving fig to "+'plots/local_dynamical_stability_NR'+str(int(NR))+'_NS'+str(int(NS))+'_exp_fit_fixed_connectance_'+label[k]+'.pdf')
#     fig.savefig('plots/local_dynamical_stability_NR'+str(int(NR))+'_NS'+str(int(NS))+'_exp_fit_fixed_connectance_'+label[k]+'.pdf')
#     plt.close(k)
