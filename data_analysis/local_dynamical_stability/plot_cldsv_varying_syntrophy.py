import numpy as np
import matplotlib.pyplot as plt
import consumer_resource_data_analysis as cf
from consumer_resource_data_analysis import alpha_mode, label, alpha0, all_nestedness, all_connectance, alpha_mode_colours, nestedness_label, connectance_label, nest_colours, conn_colours
import os
import matplotlib.tri as tr
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit

import copy
#
# filename = 'local_dynamical_stability/local_dynamical_stability_NR50_NS25_full_rank_opt_consumption_mat_NR50_NS25'
# optimal_LRI_folder='optimal_LRI_Nr50_Nc25'
# consumption_matrix_folder='optimal_matrices/consumption/Nr50_Nc25'
# matrix_set='S_{50}'

filename = 'local_dynamical_stability/all_mat_local_dynamical_stability_NR25_NS25_100_points_full_rank_opt_consumption_mat_NR25_NS25'
optimal_LRI_folder='optimal_LRI_Nr25_Nc25'
consumption_matrix_folder='optimal_matrices/consumption/Nr25_Nc25'
matrix_set='S_{25}'

# filename='local_dynamical_stability/other_LRI_local_dynamical_stability_NR25_NS25_full_rank_opt_consumption_mat_NR25_NS25'
# optimal_LRI_folder='optimal_LRI_corrected_NR25_NS25'
# consumption_matrix_folder='optimal_matrices/consumption/Nr25_Nc25'
# matrix_set='S_{25}'


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
fig.suptitle(r'$\mathcal{D}_1^{'+matrix_set+'}$ for $N_R='+str(NR)+'$, $N_S='+str(NS)+'$')
fig.savefig('plots/common_dynamical_stability_region_NR'+str(NR)+'_NS'+str(NS)+'.pdf')
#plt.show()

print("Passing to plots of local dynamical stability for each matrix ")
# now we produce a colour plot for the local dynamical stability region for all the matrices listed here
for j in range(len(local_dynamical_stability_region[0,0])):
    data = local_dynamical_stability_region[:,:,j]
    NR=data[0,0,0]
    NS=data[0,0,1]
    nestedness=data[0,0,2]
    connectance=data[0,0,3]

    fig, axs, im, levels =cf.plot_levels(data, colors, label)
    save_name='NR'+str(int(NR))+'_NS'+str(int(NS))+'_Nest'+str(nestedness)+'_Conn'+str(connectance)
    # first save figure without color bar or title
    fig.savefig('plots/local_dynamical_stability_wt_wc_region_'+save_name+'.pdf')

    # add color bar to plot
    cbar = cf.add_colorbar_to_plot_levels(fig, im, levels, alpha0)
    cbar.set_label(r'$\alpha_0$')
    fig.savefig('plots/local_dynamical_stability_wt_region_'+save_name+'.pdf')

    # add suptitle to plot
    fig.suptitle(r'Fully dynamically stable region $\mathcal{D}^{G,A}_{L,1}\left(\alpha_0\right)$ for $N_R='+str(int(NR))+', N_S='+str(int(NS))+', \kappa='+str(round(connectance,2))+', \eta='+str(nestedness)+'$')
    fig.savefig('plots/local_dynamical_stability_region_'+save_name+'.pdf')

    plt.close()

decline=[]
critical_alpha0 =[]
err_decline=[]
for j in range(len(local_dynamical_stability_region[0][0])):
    fig, ax =  plt.subplots(1, 1)
    local_critical_alpha0=[]
    local_decline=[]
    for i in range(len(local_dynamical_stability_region)):
        data=local_dynamical_stability_region[i,:,j]
        NR=data[0,0]
        NS=data[0,1]
        nestedness=cf.closest_element_in_list(data[0,2], all_nestedness)
        connectance=cf.closest_element_in_list(data[0,3], all_connectance)
        lds_volume = cf.shrink_volume_for_one_matrix(data)
        fit_function=cf.exponential_function
        take_zero_points=False
        fitted_alpha0, fitted_dyn_volume, estimated_alpha_crit, estimated_vol_zero_syntrophy, (est_decay_rate, err_decay_rate) = cf.fit_shrinkage_curve(alpha0, lds_volume, fit_function, take_zero_points)
        print(err_decay_rate)
        ax.plot(fitted_alpha0, fitted_dyn_volume,color=alpha_mode_colours[i], marker='None', linestyle='solid')
        ax.plot(alpha0, lds_volume, label=label[i], color=alpha_mode_colours[i])
        #ax.set_yscale('linear')
        local_critical_alpha0.append(estimated_alpha_crit)
        local_decline.append(est_decay_rate)
        err_decline.append(err_decay_rate)
        #print(popt)
        #ax.plot(alpha0, fitted_vol, linestyle='-', marker='None')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\alpha_0$')
    ax.set_ylabel(r'Vol$\left(\mathcal{D}^{G,A}_{L,1}(\alpha_0)\right)$')
    ax.legend()
    title=r'$N_R='+str(int(NR))+', N_S='+str(int(NS))+', \kappa='+str(round(connectance,2))+', \eta='+str(round(nestedness,2))+'$'
    ax.set_title(title)
    save_name='NR'+str(int(NR))+'_NS'+str(int(NS))+'_Nest'+str(nestedness)+'_Conn'+str(connectance)
    fig.tight_layout()
    print("Saving fig to "+'plots/size_local_dynamical_stability_region_'+save_name+'.pdf')
    fig.savefig('plots/size_local_dynamical_stability_region_'+save_name+'.pdf')
    critical_alpha0.append(local_critical_alpha0)
    decline.append(local_decline)
    plt.close()
    #plt.show()
critical_alpha0=np.array(critical_alpha0)
decline=np.array(decline)
print("Average decay rate feasibility NR=", NR, "is", np.mean(err_decline))
#
# plot a slope curve vs nestedness and connectance
critical_alpha0 = np.transpose(critical_alpha0)
decline=np.transpose(decline)

connectance = local_dynamical_stability_region[0,0][:,3]
nestedness = local_dynamical_stability_region[0,0][:,2]
NR = local_dynamical_stability_region[0,0][0,0]
NS = local_dynamical_stability_region[0,0][0,1]

for k in range(len(alpha_mode)):
    to_save = [[NR, NS, nestedness[i], connectance[i], decline[k][i]] for i in range(len(decline[k]))]
    np.savetxt(filename+'_decay_rate_'+alpha_mode[k]+'.out', np.array(to_save))


print("Passing to plots of fit")

for k in range(len(local_dynamical_stability_region)):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for j in range(len(all_nestedness)):
        nest=all_nestedness[j]
        indices = [i for i in range(len(nestedness)) if cf.closest_element_in_list(nestedness[i], all_nestedness)==nest]
        ax.plot(connectance[indices],critical_alpha0[k][indices], label=r'$\eta_G\approx'+str(nest)+'$', color=nest_colours[j])
    ax.set_xlabel(connectance_label)
    ax.set_ylabel(r'$\alpha_C^D(G,A)$')
    ax.set_title(label[k])
    ax.legend(bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout()
    print("Saving fig to "+'plots/local_dynamical_stability_NR'+str(int(NR))+'_NS'+str(int(NS))+'_critical_dynamical_syntrophy_fixed_nestedness_'+label[k]+'.pdf')
    fig.savefig('plots/local_dynamical_stability_NR'+str(int(NR))+'_NS'+str(int(NS))+'_critical_dynamical_syntrophy_fixed_nestedness_'+alpha_mode[k]+'.pdf')
    plt.close()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for j in range(len(all_connectance)):
        conn=all_connectance[j]
        indices = [i for i in range(len(connectance)) if cf.closest_element_in_list(connectance[i], all_connectance)==conn]
        sorted_indices=[ indices[a] for a in np.argsort(nestedness[indices])]
        ax.plot(nestedness[sorted_indices],critical_alpha0[k][sorted_indices], label=r'$\kappa_G\approx'+str(conn)+'$')
    ax.set_xlabel(nestedness_label)
    ax.set_ylabel(r'$\alpha_C^D(G,A)$')

    ax.set_title(label[k])
    ax.legend(bbox_to_anchor=(1.0, 1.0))

    fig.tight_layout()
    print("Saving fig to "+'plots/local_dynamical_stability_NR'+str(int(NR))+'_NS'+str(int(NS))+'_critical_dynamical_syntrophy_fixed_connectance_'+label[k]+'.pdf')
    fig.savefig('plots/local_dynamical_stability_NR'+str(int(NR))+'_NS'+str(int(NS))+'_critical_dynamical_syntrophy_fixed_connectance_'+alpha_mode[k]+'.pdf')
    plt.close()

for k in range(len(local_dynamical_stability_region)):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for j in range(len(all_nestedness)):
        nest=all_nestedness[j]
        indices = [i for i in range(len(nestedness)) if cf.closest_element_in_list(nestedness[i], all_nestedness)==nest]
        sorted_indices=[indices[a] for a in np.argsort(connectance[indices])]
        ax.plot(connectance[sorted_indices],decline[k][sorted_indices], label=r'$\eta_G\approx'+str(nest)+'$',color=nest_colours[j])
    ax.set_xlabel(connectance_label)
    ax.set_ylabel(r'$d_D(G,A)$')
    ax.set_title(label[k])
    ax.legend(bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout()
    print("Saving fig to "+'plots/local_dynamical_stability_NR'+str(int(NR))+'_NS'+str(int(NS))+'_dyn_stability_decay_rate_fixed_nestedness_'+label[k]+'.pdf')
    fig.savefig('plots/local_dynamical_stability_NR'+str(int(NR))+'_NS'+str(int(NS))+'_dyn_stability_decay_rate_fixed_nestedness_'+alpha_mode[k]+'.pdf')
    plt.close()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for j in range(len(all_connectance)):
        conn=all_connectance[j]
        indices = [i for i in range(len(connectance)) if cf.closest_element_in_list(connectance[i], all_connectance)==conn]
        sorted_indices=[ indices[a] for a in np.argsort(nestedness[indices])]
        ax.plot(nestedness[sorted_indices],decline[k][sorted_indices], label=r'$\kappa_G\approx'+str(conn)+'$',color=nest_colours[j])
    ax.set_xlabel(nestedness_label)
    ax.set_ylabel(r'$d_D(G,A)$')

    ax.set_title(label[k])
    ax.legend(bbox_to_anchor=(1.0, 1.0))

    fig.tight_layout()
    print("Saving fig to "+'plots/local_dynamical_stability_NR'+str(int(NR))+'_NS'+str(int(NS))+'_dyn_stability_decay_rate_fixed_connectance_'+label[k]+'.pdf')
    fig.savefig('plots/local_dynamical_stability_NR'+str(int(NR))+'_NS'+str(int(NS))+'_dyn_stability_decay_rate_fixed_connectance_'+alpha_mode[k]+'.pdf')
    plt.close()

# same plots but deviation from FC case
for k in range(len(local_dynamical_stability_region)):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for j in range(len(all_nestedness)):
        nest=all_nestedness[j]
        indices = [i for i in range(len(nestedness)) if cf.closest_element_in_list(nestedness[i], all_nestedness)==nest]
        sorted_indices=[indices[a] for a in np.argsort(connectance[indices])]
        ax.plot(connectance[sorted_indices],1-decline[k][sorted_indices]/decline[0][sorted_indices], label=r'$\eta_G\approx'+str(nest)+'$',color=nest_colours[j])
    ax.set_xlabel(connectance_label)
    ax.set_ylabel(r'$1-d_D(G,A)/d_D(G,$FC$)$')
    ax.set_title(label[k])
    ax.legend(bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout()
    print("Saving fig to "+'plots/local_dynamical_stability_NR'+str(int(NR))+'_NS'+str(int(NS))+'_dyn_stability_decay_rate_away_from_FC_fixed_nestedness_'+label[k]+'.pdf')
    fig.savefig('plots/local_dynamical_stability_NR'+str(int(NR))+'_NS'+str(int(NS))+'_dyn_stability_decay_rate_away_from_FC_fixed_nestedness_'+alpha_mode[k]+'.pdf')
    plt.close()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for j in range(len(all_connectance)):
        conn=all_connectance[j]
        indices = [i for i in range(len(connectance)) if cf.closest_element_in_list(connectance[i], all_connectance)==conn]
        sorted_indices=[ indices[a] for a in np.argsort(nestedness[indices])]
        ax.plot(nestedness[sorted_indices],1-decline[k][sorted_indices]/decline[0][sorted_indices], label=r'$\kappa_G\approx'+str(conn)+'$', color=conn_colours[j])
    ax.set_xlabel(nestedness_label)
    ax.set_ylabel(r'$1-d_D(G,A)/d_D(G,$FC$)$')

    ax.set_title(label[k])
    ax.legend(bbox_to_anchor=(1.0, 1.0))

    fig.tight_layout()
    print("Saving fig to "+'plots/local_dynamical_stability_NR'+str(int(NR))+'_NS'+str(int(NS))+'_dyn_stability_decay_rate_away_from_FC_fixed_connectance_'+label[k]+'.pdf')
    fig.savefig('plots/local_dynamical_stability_NR'+str(int(NR))+'_NS'+str(int(NS))+'_dyn_stability_decay_rate_away_from_FC_syntrophy_fixed_connectance_'+alpha_mode[k]+'.pdf')
    plt.close()
