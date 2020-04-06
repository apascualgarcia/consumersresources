import numpy as np
import matplotlib.pyplot as plt
import common_functions as cf
import os
import matplotlib.tri as tr
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit

import copy

alpha_mode=['random_structure', 'no_release_when_eat', 'optimal_matrix']
label=['fully connected', 'intraspecific syntrophy restricted', 'optimal LRI']
filename = 'feasibility/feasibility_NR50_NS25_full_rank_opt_consumption_mat_NR50_NS25'
alpha0=[0, 1.3e-3, 2.6e-3, 3.9e-3, 5.2e-3, 6.5e-3, 7.8e-3, 9.1e-3, 1.04e-2, 1.4e-2]
all_nestedness=[0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]
all_connectance=[0.08, 0.13, 0.18, 0.23, 0.28, 0.33, 0.38, 0.43]
cmap = plt.cm.get_cmap('jet_r')
colors = [cmap(i/10) for i in range(len(alpha0))]
alpha0_mode_colours=['blue', 'orange', 'green']

# feasibility region[alpha_mode][alpha0][connectance][nestedness][gamma0][S0] contains the feasibility of said point
feasibility_region = []
for al_mo in alpha_mode:
    local_vector=[]
    for a in alpha0:
        #file = filename+'_'+al_mo+'_optimal_LRI_alpha0='+str(a)
        #cf.remove_strings_from_file('optimal_matrices/consumption/Nr25_Nc25', file)
        file = filename+'_'+al_mo+'_optimal_LRI_Nr50_Nc25_alpha0='+str(a)+'_filtered.out'
        #print('Loading ', file)
        local_data=np.loadtxt(file)
        local_vector.append(local_data)
    feasibility_region.append(local_vector)
feasibility_region = np.array(feasibility_region)
alpha0=np.array(alpha0)


# Plot common feasibility volume
for k in range(len(feasibility_region)):
    fig=plt.figure(k)
    ax=fig.add_subplot(111)
    ff_indices=[]
    data=feasibility_region[0][0]
    gamma0=data[:, 4::3]
    S0=data[:, 5::3]
    Npoints=len(S0[0])
    for l in range(len(feasibility_region[k])):
        data=feasibility_region[k,l]
        NR=data[:,0]
        NS=data[:,1]
        nestedness=data[:,2]
        connectance=data[:,3]
        feasability=data[:, 6::3]
        # find the fully feasible_indices for this alpha0 and this alpha mode
        full_feas_indices=[j for j in range(len(feasability[0])) if feasability[0,j]==1.]
        for m in range(1,len(feasability)):
            oldffi=copy.deepcopy(full_feas_indices)
            full_feas_indices=[n for n in range(len(feasability[m])) if (feasability[m,n]==1. and n in oldffi)]
        in_feas_vol=[]
        for m in range(len(feasability[0])):
            if(m in full_feas_indices):
                in_feas_vol.append(1.)
            else:
                in_feas_vol.append(0.)
        ff_indices.append(full_feas_indices)

    feasibility_level=[]
    for i in range(Npoints):
        level=0
        j=0
        exit=False
        while(not(exit) and j < len(feasibility_region[k])):
            if i in ff_indices[j]:
                level+=1
            else:
                exit=True
            j+=1
        feasibility_level.append(level)
    max_level=max(feasibility_level)
    levels=[i for i in range(max_level+1)]

    isbad=np.less(feasibility_level, 0.5)
    triang = tr.Triangulation(gamma0[0], S0[0])
    #ax.tricontour(triang, feasibility_level, levels=levels, colors='black')
    mask = np.all(np.where(isbad[triang.triangles], True, False), axis=1)
    triang.set_mask(mask)
    im=ax.tricontourf(triang, feasibility_level, levels=levels, colors=colors)

    ax.set_xlim(0., 1.)
    ax.set_ylim(0., 1.)
    ax.set_xlabel(r'$\gamma_0$')
    ax.set_ylabel(r'$S_0$')
    cbar=fig.colorbar(im)
    cbar.set_ticks([a +0.5 for a in levels[0:max_level+2]])
    cbar.set_ticklabels(alpha0)
    cbar.set_label(r'$\alpha_0$')
    ax.set_title(r'$\mathcal{F}^*_1$ for $N_R='+str(int(NR[0]))+', N_S='+str(int(NS[0]))+'$ '+label[k])
    fig.tight_layout()
    #fig.savefig('plots/common_feasibility_volume_NR'+str(int(NR[0]))+'_NS'+str(int(NS[0]))+'_varying_syntrophy_'+alpha_mode[k]+'.pdf')
    plt.close(k)
#
# now we produce a similar plot for all the matrices listed here
for j in range(len(feasibility_region[0,0])):
#for j in range(1):

    fig, axs = plt.subplots(1, 3, sharey=True, sharex=True, figsize=(10.5,4.5))
    feasible_levels=[]

    data = feasibility_region[0,:,0]
    gamma0=data[:, 4::3]
    S0=data[:, 5::3]
    NR=data[0,0]
    NS=data[0,1]

    nestedness=data[0,2]
    connectance=data[0,3]

    # first get the different levels for the three alpha modes
    for i in range(len(feasibility_region)):
        nestedness=data[0,2]
        connectance=data[0,3]
        data = feasibility_region[i,:,j]
        feasible_levels.append(cf.data_levels(data))
    feasible_levels=np.array(feasible_levels)
    max_level=np.amax(feasible_levels)
    # the unfeasible level is left blank
    levels=[i for i in range(max_level+2)]

    # then we found the levels we can plot the data
    triang = tr.Triangulation(gamma0[0], S0[0])
    for i in range(len(feasibility_region)):
        to_plot = feasible_levels[i]
        axs[i].set_aspect('equal')
        im = axs[i].tricontourf(triang, to_plot, levels=levels, colors=colors)
        axs[i].set_xlabel(r'$\gamma_0$')
        axs[i].set_xticks([0, 0.5, 1])
        axs[i].set_xticklabels([0, 0.5, 1])
        axs[i].set_title(label[i])

    axs[0].set_ylabel(r'$S_0$')
    axs[0].set_yticks([0, 0.5, 1])
    axs[0].set_yticklabels([0, 0.5, 1])
    fig.subplots_adjust(bottom=0.2, top=0.95)
    save_name='NR'+str(int(NR))+'_NS'+str(int(NS))+'_Nest'+str(nestedness)+'_Conn'+str(connectance)
    #fig.savefig('plots/feasibility_region_wt_wc_'+save_name+'.pdf')
    cbar_ax = fig.add_axes([0.125, 0.15, 0.75, 0.02])
    cbar=fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.set_label(r'$\alpha_0$')
    cbar.set_ticks([a +0.5 for a in levels])
    cbar.set_ticklabels(alpha0)
    #fig.savefig('plots/feasibility_region_wt_'+save_name+'.pdf')
    fig.suptitle(r'Fully feasible region $\mathcal{F}^G_1$ for $N_R='+str(int(NR))+', N_S='+str(int(NS))+', \kappa='+str(round(connectance,2))+', \eta='+str(nestedness)+'$')
    #fig.savefig('plots/feasibility_region_'+save_name+'.pdf')
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
        ax.plot(alpha0, fitted_vol, marker='None', linestyle='solid', linewidth=2, color=alpha0_mode_colours[i])
        ax.plot(alpha0, feas_volume, label=label[i], markersize=10, linestyle='None', markeredgewidth=3, color=alpha0_mode_colours[i])
    ax.set_yscale('log')
    ax.set_xlabel(r'$\alpha_0$')
    ax.set_ylabel(r'Vol$(\mathcal{F}^G_1(\alpha_0))$ (normalized)')
    ax.set_ylim(-0.02,0.55)
    ax.legend()
    title=r'$N_R='+str(int(NR))+', N_S='+str(int(NS))+', \kappa='+str(round(connectance,2))+', \eta='+str(round(nestedness,2))+'$'
    ax.set_title(title)
    save_name='NR'+str(int(NR))+'_NS'+str(int(NS))+'_Nest'+str(nestedness)+'_Conn'+str(connectance)
    fig.tight_layout()
    #fig.savefig('plots/size_feasibility_region_'+save_name+'.pdf')
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

    #fig.savefig('plots/feasibility_NR'+str(int(NR))+'_NS'+str(int(NS))+'_feasibility_decay_rate_fixed_nestedness_'+alpha_mode[k]+'.pdf')
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
    #fig.savefig('plots/feasibility_NR'+str(int(NR))+'_NS'+str(int(NS))+'_feasibility_decay_rate_fixed_connectance_'+alpha_mode[k]+'.pdf')
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
    #fig.savefig('plots/feasibility_NR'+str(int(NR))+'_NS'+str(int(NS))+'_critical_alpha0_fixed_nestedness_'+alpha_mode[k]+'.pdf')
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
    #fig.savefig('plots/feasibility_NR'+str(int(NR))+'_NS'+str(int(NS))+'_critical_alpha0_fixed_connectance_'+alpha_mode[k]+'.pdf')
    plt.close()
