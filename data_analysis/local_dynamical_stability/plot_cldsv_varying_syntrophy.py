import numpy as np
import matplotlib.pyplot as plt
import common_functions as cf
import os
import matplotlib.tri as tr
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit

import copy

alpha_mode=['random_structure', 'no_release_when_eat', 'optimal_matrix']
label=['fully connected', 'no release when eat', 'optimal LRI']
filename = 'local_dynamical_stability/local_dynamical_stability_NR25_NS25_full_rank_opt_consumption_mat_NR25_NS25'
alpha0=[0, 1.3e-3, 2.6e-3, 3.9e-3, 5.2e-3, 6.5e-3, 7.8e-3, 9.1e-3, 1.04e-2, 1.4e-2]
all_nestedness=[0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]
all_connectance=[0.08, 0.13, 0.18, 0.23, 0.28, 0.33, 0.38, 0.43]

cmap = plt.cm.get_cmap('jet_r')
colors = [cmap(i/10) for i in range(len(alpha0))]

# local_dynamical_stability region[alpha_mode][alpha0][connectance][nestedness][gamma0][S0] contains the local_dynamical_stability of said point
local_dynamical_stability_region = []
for al_mo in alpha_mode:
    local_vector=[]
    for a in alpha0:
        file = filename+'_'+al_mo+'_optimal_LRI_alpha0='+str(a)+'_filtered.out'
        local_data=np.loadtxt(file)
        local_vector.append(local_data)
    local_dynamical_stability_region.append(local_vector)
local_dynamical_stability_region = np.array(local_dynamical_stability_region)

alpha0=np.array(alpha0)


# Plot common local_dynamical_stability volume
for k in range(len(local_dynamical_stability_region)):
    fig=plt.figure(k)
    ax=fig.add_subplot(111)
    ff_indices=[]
    data=local_dynamical_stability_region[0][0]
    gamma0=data[:, 4::3]
    S0=data[:, 5::3]
    Npoints=len(S0[0])
    for l in range(len(local_dynamical_stability_region[k])):
        data=local_dynamical_stability_region[k][l]
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

    local_dynamical_stability_level=[]
    for i in range(Npoints):
        level=0
        j=0
        exit=False
        while(not(exit) and j < len(local_dynamical_stability_region[k])):
            if i in ff_indices[j]:
                level+=1
            else:
                exit=True
            j+=1
        local_dynamical_stability_level.append(level)
    max_level=max(local_dynamical_stability_level)
    levels=[i for i in range(max_level+1)]

    isbad=np.less(local_dynamical_stability_level, 0.5)
    triang = tr.Triangulation(gamma0[0], S0[0])
    #ax.tricontour(triang, local_dynamical_stability_level, levels=levels, colors='black')
    mask = np.all(np.where(isbad[triang.triangles], True, False), axis=1)
    triang.set_mask(mask)
    im=ax.tricontourf(triang, local_dynamical_stability_level, levels=levels, colors=colors)

    ax.set_xlim(0., 1.)
    ax.set_ylim(0., 1.)
    ax.set_xlabel(r'$\gamma_0$')
    ax.set_ylabel(r'$S_0$')
    cbar=fig.colorbar(im)
    cbar.set_ticks([a +0.5 for a in levels[0:max_level+2]])
    cbar.set_ticklabels(alpha0)
    cbar.set_label(r'$\alpha_0$')
    ax.set_title(r'$\mathcal{D}_{L,1}^*$ for $N_R=25, N_S=25$ '+label[k])
    fig.tight_layout()
    fig.savefig("plots/common_local_dynamical_stability_volume_varying_syntrophy_"+alpha_mode[k]+'.pdf')
    plt.close(k)
# now we produce a similar plot for all the matrices listed here
for j in range(len(local_dynamical_stability_region[0,0])):
#for j in range(1):

    fig, axs = plt.subplots(1, 3, sharey=True, sharex=True, figsize=(10.5,4.5))
    lds_levels=[]

    data = local_dynamical_stability_region[0,:,0]
    gamma0=data[:, 4::3]
    S0=data[:, 5::3]
    NR=data[0,0]
    NS=data[0,1]

    nestedness=data[0,2]
    connectance=data[0,3]

    # first get the different levels for the three alpha modes
    for i in range(len(local_dynamical_stability_region)):
        data = local_dynamical_stability_region[i,:,j]
        nestedness=data[0,2]
        connectance=data[0,3]
        lds_levels.append(cf.data_levels(data))
    lds_levels=np.array(lds_levels)
    max_level=np.amax(lds_levels)
    # the unfeasible level is left blank
    levels=[i for i in range(max_level+2)]

    # then we found the levels we can plot the data
    triang = tr.Triangulation(gamma0[0], S0[0])
    for i in range(len(local_dynamical_stability_region)):
        to_plot = lds_levels[i]
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
    cbar_ax = fig.add_axes([0.125, 0.15, 0.75, 0.02])
    cbar=fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.set_label(r'$\alpha_0$')
    cbar.set_ticks([a +0.5 for a in levels])
    cbar.set_ticklabels(alpha0)

    fig.suptitle(r'Fully dynamically stable region $\mathcal{D}^G_{L,1}$ for $N_R='+str(int(NR))+', N_S='+str(int(NS))+', \kappa='+str(round(connectance,2))+', \eta='+str(nestedness)+'$')
    save_name='NR'+str(int(NR))+'_NS'+str(int(NS))+'_Nest'+str(nestedness)+'_Conn'+str(connectance)
    fig.savefig('plots/local_dynamical_stability_region_'+save_name+'.pdf')
    plt.close()

decline =[]
for j in range(len(local_dynamical_stability_region[0][0])):
    #fig, ax =  plt.subplots(1, 1)
    local_decline=[]
    for i in range(len(local_dynamical_stability_region)):
        data=local_dynamical_stability_region[i,:,j]
        NR=data[0,0]
        NS=data[0,1]
        nestedness=data[0,2]
        connectance=data[0,3]
        feas_volume = cf.shrink_volume_for_one_matrix(data)
        # ax.plot(alpha0, feas_volume, label=label[i])
        popt, pcov=curve_fit(cf.exp_fit, alpha0, feas_volume, maxfev=90000000)
        # fitted_vol=[cf.exp_fit(x,*popt) for x in alpha0]
        local_decline.append(popt[1])
        #print(popt)
        #ax.plot(alpha0, fitted_vol, linestyle='-', marker='None')
    # ax.set_xlabel(r'$\alpha_0$')
    # ax.set_ylabel(r'Vol$(\mathcal{V}^*(\alpha_0))$ (normalized)')
    # ax.set_ylim(-0.02,0.55)
    # ax.legend()
    # title=r'$N_R='+str(int(NR))+', N_S='+str(int(NS))+', \kappa='+str(round(connectance,2))+', \eta='+str(round(nestedness,2))+'$'
    # ax.set_title(title)
    # save_name='NR'+str(int(NR))+'_NS'+str(int(NS))+'_Nest'+str(nestedness)+'_Conn'+str(connectance)
    # fig.tight_layout()
    # fig.savefig('plots/size_local_dynamical_stability_region_'+save_name+'.pdf')
    decline.append(local_decline)
    #plt.show()
decline=np.array(decline)
#
# plot a slope curve vs nestedness and connectance
decline = np.transpose(decline)

connectance = local_dynamical_stability_region[0,0][:,3]
nestedness = local_dynamical_stability_region[0,0][:,2]

for k in range(len(local_dynamical_stability_region)):
    fig = plt.figure(k)
    ax = fig.add_subplot(111)
    for nest in all_nestedness:
        indices = [i for i in range(len(nestedness)) if cf.closest_element_in_list(nestedness[i], all_nestedness)==nest]
        ax.plot(connectance[indices],decline[k][indices], label=r'$\eta\approx'+str(nest)+'$')
    ax.set_xlabel(r'Connectance $\kappa$')
    ax.set_ylabel(r'Exponential coefficient')
    ax.set_title(label[k])
    ax.legend(bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout()
    fig.savefig('plots/local_dynamical_stability_NR25_NS25_exp_fit_fixed_nestedness_'+label[k]+'.pdf')
    plt.close(k)

    fig = plt.figure(k)
    ax = fig.add_subplot(111)
    for conn in all_connectance:
        indices = [i for i in range(len(connectance)) if cf.closest_element_in_list(connectance[i], all_connectance)==conn]
        sorted_indices=[ indices[a] for a in np.argsort(nestedness[indices])]
        ax.plot(nestedness[sorted_indices],decline[k][sorted_indices], label=r'$\kappa\approx'+str(conn)+'$')
    ax.set_xlabel(r'Nestedness $\eta$')
    ax.set_ylabel(r'Exponential coefficient')
    ax.set_title(label[k])
    ax.legend(bbox_to_anchor=(1.0, 1.0))

    fig.tight_layout()
    fig.savefig('plots/local_dynamical_stability_NR25_NS25_exp_fit_fixed_connectance_'+label[k]+'.pdf')
    plt.close(k)
