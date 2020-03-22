import numpy as np
import matplotlib.pyplot as plt
from common_functions import plot_points_as_area, get_fit_cfr
import os
import matplotlib.tri as tr

import copy

alpha_mode=['random_structure', 'no_release_when_eat', 'optimal_matrix']
label=['random structure', 'no release when eat', 'optimal LRI']
filename = 'data_output/feasibility_NR25_NS25_full_rank_opt_consumption_mat_NR25_NS25'
alpha0=[0, 0.0013, 0.0026, 0.0039, 0.0052, 0.0065, 0.0078, 0.0091, 0.0104, 0.014]
# feasibility region[alpha_mode][alpha0][connectance][nestedness][gamma0][S0] contains the feasibility of said point
feasibility_region = []
for al_mo in alpha_mode:
    local_vector=[]
    for a in alpha0:
        file = filename+'_'+al_mo+'_optimal_LRI_alpha0='+str(a)+'_filtered.out'
        local_data=np.loadtxt(file)
        local_vector.append(local_data)
    feasibility_region.append(local_vector)


for k in range(len(feasibility_region)):
    fig=plt.figure(k)
    ax=fig.add_subplot(111)
    ff_indices=[]
    data=feasibility_region[0][0]
    gamma0=data[:, 4::3]
    S0=data[:, 5::3]
    Npoints=len(S0[0])
    for l in range(len(feasibility_region[k])):
        data=feasibility_region[k][l]
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
    im=ax.tricontourf(triang, feasibility_level, levels=levels, cmap='jet_r')

    ax.set_xlim(0., 1.)
    ax.set_ylim(0., 1.)
    ax.set_xlabel(r'$\gamma_0$')
    ax.set_ylabel(r'$S_0$')
    cbar=fig.colorbar(im)
    cbar.set_ticks(levels[0:max_level])
    cbar.set_ticklabels(alpha0[0:max_level])
    ax.set_title(r'$N_R=25, N_S=25$ '+label[k])
    fig.tight_layout()
    fig.savefig("plots/common_feasibility_volume_varying_syntrophy_"+str(alpha_mode)+'.pdf')
