import numpy as np
import matplotlib.pyplot as plt
import common_functions as cf
import os
import sys
import matplotlib.tri as tr
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit

import copy

np.set_printoptions(threshold=sys.maxsize)

alpha_mode=['random_structure', 'no_release_when_eat', 'optimal_matrix']
label=['fully connected', 'no intraspecific syntrophy', 'optimal LRI']
filename_dynamical = 'local_dynamical_stability/local_dynamical_stability_NR50_NS25_full_rank_opt_consumption_mat_NR50_NS25'
filename_feasible = 'feasibility/feasibility_NR50_NS25_full_rank_opt_consumption_mat_NR50_NS25'
alpha0=[0, 1.3e-3, 2.6e-3, 3.9e-3, 5.2e-3, 6.5e-3, 7.8e-3, 9.1e-3, 1.04e-2, 1.4e-2]
all_nestedness=[0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]
all_connectance=[0.08, 0.13, 0.18, 0.23, 0.28, 0.33, 0.38, 0.43]

cmap = plt.cm.get_cmap('jet_r')
colors = [cmap(i/10) for i in range(len(alpha0))]

# first load both local dynamical stability and feasibility data
local_dynamical_stability_region = []
for al_mo in alpha_mode:
    local_vector=[]
    for a in alpha0:
        file = filename_dynamical+'_'+al_mo+'_optimal_LRI_Nr50_Nc25_alpha0='+str(a)+'_filtered.out'
        local_data=np.loadtxt(file)
        local_vector.append(local_data)
    local_dynamical_stability_region.append(local_vector)
local_dynamical_stability_region = np.array(local_dynamical_stability_region)

feasibility_region = []
for al_mo in alpha_mode:
    local_vector=[]
    for a in alpha0:
        file = filename_feasible+'_'+al_mo+'_optimal_LRI_Nr50_Nc25_alpha0='+str(a)+'_filtered.out'
        local_data=np.loadtxt(file)
        local_vector.append(local_data)
    feasibility_region.append(local_vector)
feasibility_region = np.array(feasibility_region)
alpha0=np.array(alpha0)

# for all matrices, at each mode and each alpha0, find points which are not both
# in the feasibility and local dynamical stability regions
for k in range(len(feasibility_region[0,0])):
#for k in range(1):
    fig, axs = plt.subplots(1, 3, sharey=True, sharex=True, figsize=(10.5,4.5))
    for i in range(len(feasibility_region)):
    #for i in range(1):
        ratio_lds_feasible=[]
        for j in range(len(feasibility_region[i])):
        #for j in range(1):
            feasibility_data=feasibility_region[i,j,k]

            NR = feasibility_data[0]
            NS = feasibility_data[1]
            nestedness=feasibility_data[2]
            connectance=feasibility_data[3]

            gamma0 = feasibility_data[4::3]
            S0 = feasibility_data[5::3]
            feasibility = feasibility_data[6::3]

            lds_data=local_dynamical_stability_region[i,j,k]
            lds = lds_data[6::3]

            # find the range of indices which are fully feasible and locally dynamically stable
            feasible_indices=[]
            lds_indices=[]
            for l in range(len(feasibility)):
                if feasibility[l]==1.:
                    feasible_indices.append(l)
                if lds[l]==1.:
                    lds_indices.append(l)
            lds_indices=[i for i in lds_indices if i in feasible_indices]
            only_feasible = [l for l in feasible_indices if l not in lds_indices]
            if len(feasible_indices)>0:
                ratio = len(lds_indices)/len(feasible_indices)
                ratio_lds_feasible.append(ratio)
        axs[i].set_xlabel(r'$\alpha_0$')
        axs[i].set_yscale('linear')
        axs[i].plot(alpha0[0:len(ratio_lds_feasible)], ratio_lds_feasible,markersize=10, linewidth=2.5, markeredgewidth=3)
        axs[i].set_title(label[i])
        axs[i].set_xlim(-0.002, max(alpha0)*1.1)
    axs[0].set_ylabel(r'Vol$\left(\mathcal{D}_{L,1}^{G,A}(\alpha_0)\right)$/Vol$\left(\mathcal{F}_{1}^{G,A}(\alpha_0)\right)$')
    title=r'$N_R='+str(int(NR))+', N_S='+str(int(NS))+', \kappa='+str(round(connectance,2))+', \eta='+str(round(nestedness,2))+'$'
    save_name='NR'+str(int(NR))+'_NS'+str(int(NS))+'_Nest'+str(nestedness)+'_Conn'+str(connectance)
    fig.subplots_adjust(top=0.7)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig('plots/feasibility_vs_lds_'+save_name+'.pdf')
    plt.close()

# do the same computation for the common feasible and dynamical volumes
# Plot common feasibility volume
for k in range(len(feasibility_region)):
    fig=plt.figure(k)
    ax=fig.add_subplot(111)

    ## ff_indices[l] contains all fully feasible indices at alpha0[l]
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

    ## fds_indices[l] contains all fully dynamically stable indices at alpha0[l]
    fds_indices=[]
    data=local_dynamical_stability_region[0][0]
    gamma0=data[:, 4::3]
    S0=data[:, 5::3]
    Npoints=len(S0[0])
    for l in range(len(local_dynamical_stability_region[k])):
        data=local_dynamical_stability_region[k,l]
        NR=data[:,0]
        NS=data[:,1]
        nestedness=data[:,2]
        connectance=data[:,3]
        local_dynamical_stability=data[:, 6::3]
        # find the fully lds indices for this alpha0 and this alpha mode
        full_lds_indices=[j for j in range(len(local_dynamical_stability[0])) if local_dynamical_stability[0,j]==1.]
        for m in range(1,len(feasability)):
            oldfldsi=copy.deepcopy(full_lds_indices)
            full_lds_indices=[n for n in range(len(local_dynamical_stability[m])) if (local_dynamical_stability[m,n]==1. and n in oldfldsi)]
        fds_indices.append(full_lds_indices)

    for l in range(len(local_dynamical_stability_region[k])):
        fds_indices[l]=[i for i in fds_indices[l] if i in ff_indices[l]]
        only_feasible=[i for i in ff_indices[l] if i not in fds_indices[l]]
        data=local_dynamical_stability_region[k,l]
        print("For alpha0 = ", alpha0[l], " only feasible are ", only_feasible)
