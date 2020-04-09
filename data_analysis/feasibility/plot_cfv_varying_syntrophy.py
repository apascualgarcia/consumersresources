import numpy as np
import matplotlib.pyplot as plt
import consumer_resource_data_analysis as cf
from consumer_resource_data_analysis import alpha_mode, label, alpha0, all_nestedness, all_connectance,alpha_mode_colours
import os
import matplotlib.tri as tr
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit

import copy

filename = 'feasibility/feasibility_NR25_NS25_100_points_full_rank_opt_consumption_mat_NR25_NS25'
optimal_LRI_folder='optimal_LRI'
consumption_matrix_folder='optimal_matrices/consumption/Nr25_Nc25'

cmap = plt.cm.get_cmap('jet_r')
colors = [cmap(i/10) for i in range(len(alpha0))]

# FIRST LOAD DATA
# feasibility region[alpha_mode][alpha0][connectance][nestedness][gamma0][S0] contains the feasibility of said point
cf.filter_data(alpha_mode, alpha0, filename, optimal_LRI_folder, consumption_matrix_folder)
feasibility_region = cf.load_data_region(alpha_mode, alpha0, filename, optimal_LRI_folder)
alpha0=np.array(alpha0)

# Plot common feasibility volume
NR=int(feasibility_region[0,0,0,0])
NS=int(feasibility_region[0,0,0,1])
fig, axs, im, levels = cf.plot_common_region(feasibility_region, alpha_mode, colors, label)
cbar = cf.add_colorbar_to_plot_levels(fig, im, levels, alpha0)
cbar.set_label(r'$\alpha_0$')
fig.suptitle(r'$\mathcal{F}_1^{S_M}$ for $N_R='+str(NR)+'$, $N_S='+str(NS)+'$')
fig.savefig('plots/common_feasibility_region_NR'+str(NR)+'_NS'+str(NS)+'.pdf')

# plot feasibility levels for all matrices
for j in range(len(feasibility_region[0,0])):
    data = feasibility_region[:,:,j]
    NR = data[0,0,0]
    NS = data[0,0,1]
    nestedness = data[0,0,2]
    connectance = data[0,0,3]

    fig, axs, im, levels = cf.plot_levels(data, colors, label)
    save_name='NR'+str(int(NR))+'_NS'+str(int(NS))+'_Nest'+str(nestedness)+'_Conn'+str(connectance)
    # first save figure without color bar or title
    fig.savefig('plots/feasibility_region_wt_wc_'+save_name+'.pdf')

    # add colorbar to plot
    cbar = cf.add_colorbar_to_plot_levels(fig, im, levels, alpha0)
    cbar.set_label(r'$\alpha_0$')
    fig.savefig('plots/feasibility_region_wt_'+save_name+'.pdf')

    # add suptitle to plot
    fig.suptitle(r'Fully feasible region $\mathcal{F}^{G,A}_1$ for $N_R='+str(int(NR))+', N_S='+str(int(NS))+', \kappa='+str(round(connectance,2))+', \eta='+str(nestedness)+'$')
    fig.savefig('plots/feasibility_region_'+save_name+'.pdf')
    plt.close()


decline =[]
critical_alpha0=[]
for j in range(len(feasibility_region[0][0])):
    fig, ax =  plt.subplots(1, 1)
    local_decline=[]
    local_critical=[]
    for i in range(len(feasibility_region)):
        data=feasibility_region[i,:,j]
        NR=data[0,0]
        NS=data[0,1]
        nestedness=data[0,2]
        connectance=data[0,3]
        feas_volume = cf.shrink_volume_for_one_matrix(data)
        fitted_vol, popt, pcov=cf.fit_data(cf.exponential_function, alpha0, feas_volume)
        crit_alpha0, error=cf.zero_from_fit(cf.exponential_function, popt, pcov)
        local_decline.append(popt[1])
        local_critical.append(crit_alpha0)
        ax.plot(alpha0, fitted_vol, marker='None', linestyle='solid', linewidth=2, color=alpha_mode_colours[i])
        ax.plot(alpha0, feas_volume, label=label[i], markersize=10, linestyle='None', markeredgewidth=3, color=alpha_mode_colours[i])
    ax.set_yscale('log')
    ax.set_xlabel(r'$\alpha_0$')
    ax.set_ylabel(r'Vol$(\mathcal{F}^G_1(\alpha_0))$ (normalized)')
    ax.set_ylim(-0.02,0.55)
    ax.legend()
    title=r'$N_R='+str(int(NR))+', N_S='+str(int(NS))+', \kappa='+str(round(connectance,2))+', \eta='+str(round(nestedness,2))+'$'
    ax.set_title(title)
    save_name='NR'+str(int(NR))+'_NS'+str(int(NS))+'_Nest'+str(nestedness)+'_Conn'+str(connectance)
    fig.tight_layout()
    print('Saving fig under', 'plots/size_feasibility_region_'+save_name+'.pdf')
    fig.savefig('plots/size_feasibility_region_'+save_name+'.pdf')
    decline.append(local_decline)
    critical_alpha0.append(local_critical)
    #plt.show()
    plt.close()

decline=np.array(decline)
critical_alpha0=np.array(critical_alpha0)

decline = np.transpose(decline)
critical_alpha0=np.transpose(critical_alpha0)



# plot a feasibility decay rate curve vs nestedness and connectance
connectance = feasibility_region[0,0][:,3]
nestedness = feasibility_region[0,0][:,2]
NR = feasibility_region[0,0][0,0]
NS = feasibility_region[0,0][0,1]

for k in range(len(alpha_mode)):
    to_save = [[NR, NS,nestedness[i], connectance[i], decline[k][i]] for i in range(len(decline[k]))]
    np.savetxt(filename+'_decay_rate_'+alpha_mode[k]+'.out', np.array(to_save))



print(connectance)
print(nestedness)
print(decline)
ylim = (0, np.max(decline)*1.1)

for k in range(len(feasibility_region)):
    fig = plt.figure(k)
    ax = fig.add_subplot(111)
    for nest in all_nestedness:
        indices = [i for i in range(len(nestedness)) if cf.closest_element_in_list(nestedness[i], all_nestedness)==nest]
        ax.plot(connectance[indices],decline[k][indices], label=r'$\eta\approx'+str(nest)+'$',markersize=10, linewidth=2.5, markeredgewidth=3)
    ax.set_ylim(ylim)
    ax.set_xlabel(r'Connectance $\kappa$')
    ax.set_ylabel(r'Feasibility decay rate $d_F$')
    ax.set_title(label[k])
    ax.legend(bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout()

    fig.savefig('plots/feasibility_NR'+str(int(NR))+'_NS'+str(int(NS))+'_feasibility_decay_rate_fixed_nestedness_'+alpha_mode[k]+'.pdf')
    plt.close(k)

    fig = plt.figure(k)
    ax = fig.add_subplot(111)
    for conn in all_connectance:
        indices = [i for i in range(len(connectance)) if cf.closest_element_in_list(connectance[i], all_connectance)==conn]
        sorted_indices=[ indices[a] for a in np.argsort(nestedness[indices])]
        ax.plot(nestedness[sorted_indices],decline[k][sorted_indices], label=r'$\kappa\approx'+str(conn)+'$',markersize=10, linewidth=2.5, markeredgewidth=3)
    ax.set_ylim(ylim)
    ax.set_xlabel(r'Ecological overlap $\eta$')
    ax.set_ylabel(r'Feasibility decay rate $d_F$')
    ax.legend(bbox_to_anchor=(1.0, 1.0))
    ax.set_title(label[k])
    fig.tight_layout()
    fig.savefig('plots/feasibility_NR'+str(int(NR))+'_NS'+str(int(NS))+'_feasibility_decay_rate_fixed_connectance_'+alpha_mode[k]+'.pdf')
    plt.close(k)

ylim = (np.min(critical_alpha0)*0.9, np.max(critical_alpha0)*1.1)
# plot critical alpha0 curve vs nestedness and connectance
for k in range(len(feasibility_region)):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for nest in all_nestedness:
        indices = [i for i in range(len(nestedness)) if cf.closest_element_in_list(nestedness[i], all_nestedness)==nest]
        ax.plot(connectance[indices],critical_alpha0[k][indices], label=r'$\eta\approx'+str(nest)+'$',markersize=10, linewidth=2.5, markeredgewidth=3)
    ax.set_ylim(ylim)
    ax.set_xlabel(r'Connectance $\kappa$')
    ax.set_ylabel(r'Critical feasible syntrophy $\alpha_0^F$')
    ax.set_title(label[k])
    ax.legend(bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout()
    fig.savefig('plots/feasibility_NR'+str(int(NR))+'_NS'+str(int(NS))+'_critical_alpha0_fixed_nestedness_'+alpha_mode[k]+'.pdf')
    plt.close()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for conn in all_connectance:
        indices = [i for i in range(len(connectance)) if cf.closest_element_in_list(connectance[i], all_connectance)==conn]
        sorted_indices=[ indices[a] for a in np.argsort(nestedness[indices])]
        ax.plot(nestedness[sorted_indices],critical_alpha0[k][sorted_indices], label=r'$\kappa\approx'+str(conn)+'$',markersize=10, linewidth=2.5, markeredgewidth=3)
    ax.set_ylim(ylim)
    ax.set_xlabel(r'Ecological overlap $\eta$')
    ax.set_ylabel(r'Critical feasible syntrophy $\alpha_0^F$')
    ax.legend(bbox_to_anchor=(1.0, 1.0))
    ax.set_title(label[k])
    fig.tight_layout()
    fig.savefig('plots/feasibility_NR'+str(int(NR))+'_NS'+str(int(NS))+'_critical_alpha0_fixed_connectance_'+alpha_mode[k]+'.pdf')
    plt.close()
