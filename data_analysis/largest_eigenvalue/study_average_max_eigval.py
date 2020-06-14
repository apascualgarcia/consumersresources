import numpy as np
import matplotlib.pyplot as plt
import consumer_resource_data_analysis as cf
import os
import matplotlib.tri as tr
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit
from matplotlib.colors import LogNorm
from consumer_resource_data_analysis import alpha_mode, alpha_mode_colours,label, alpha0, all_nestedness, all_connectance, N_alphamodes
import copy

filename = 'data_shared/data_output_11-Jun_largest_eigenvalue/largest_eigenvalue_NR25_NS25_full_rank_opt_consumption_mat_NR25_NS25'
optimal_LRI_folder='optimal_LRI_Nr25_Nc25'
consumption_matrix_folder='optimal_matrices/consumption/Nr25_Nc25'

# filename = 'largest_eigenvalue/other_LRI_largest_eigenvalue_NR25_NS25_full_rank_opt_consumption_mat_NR25_NS25'
# optimal_LRI_folder='optimal_LRI_corrected_NR25_NS25'
# consumption_matrix_folder='optimal_matrices/consumption/Nr25_Nc25'
#

Npoints = 10000
square_size=8


cmap = plt.cm.get_cmap('jet_r')
colors = [cmap(i/10) for i in range(len(alpha0))]

# largest_eigenvalue region[alpha_mode][alpha0][connectance][nestedness][gamma0][S0] contains the largest_eigenvalue of said point
largest_eigenvalue_region = []
cf.filter_data(alpha_mode, alpha0, filename, optimal_LRI_folder, consumption_matrix_folder)
for al_mo in alpha_mode:
    local_vector=[]
    for a in alpha0:
        file = filename+'_'+al_mo+'_'+optimal_LRI_folder+'_alpha0='+str(a)+'_filtered.out'
        local_data=np.loadtxt(file, dtype=complex)
        print('Loading file', file, 'which contains the data of ', len(local_data), 'matrices')
        local_vector.append(local_data)
    largest_eigenvalue_region.append(local_vector)
largest_eigenvalue_region = np.array(largest_eigenvalue_region)
print('Dimensions of tableau : ', len(largest_eigenvalue_region), 'x', len(largest_eigenvalue_region[0]),'x',len(largest_eigenvalue_region[0,0]))

connectance = np.real(largest_eigenvalue_region[0,0][:,3])
nestedness = np.real(largest_eigenvalue_region[0,0][:,2])

alpha0=np.array(alpha0)

Nmatrices=len(largest_eigenvalue_region[0,0])
# now plot largest eigenvalue observed on average among all matrices
fig = plt.figure()
ax = fig.add_subplot(111)
# lists all the matrices that contribute to the average at a given alpha0
cont_matrices=[]
for i in range(len(largest_eigenvalue_region)):

    NR=np.real(largest_eigenvalue_region[0,0,0,0])
    NS=np.real(largest_eigenvalue_region[0,0,0,1])

    av_lar_eival=[]

    lcont_matrices=[]
    for j in range(len(alpha0)):
        local_cont_matrices=[]
        lar_eival_at_alpha0=[]
        for k in range(Nmatrices):
            l_eigvals=np.ma.masked_invalid(np.real(largest_eigenvalue_region[i,j,k, 6::3]))
            l_eigvals=l_eigvals[l_eigvals.mask==False]
            # if not fully unfeasible get the largest eigenvalue
            if len(l_eigvals) > 0:
                lar_eival_at_alpha0.append(np.max(l_eigvals))
                local_cont_matrices.append(k)
        if len(lar_eival_at_alpha0) > 0:
            av_lar_eival.append(np.mean(lar_eival_at_alpha0))
        else:
            av_lar_eival.append('nan')
        lcont_matrices.append(local_cont_matrices)
    av_lar_eival=np.ma.masked_invalid(av_lar_eival)
    ax.plot(alpha0, np.abs(av_lar_eival), label=label[i])
    cont_matrices.append(lcont_matrices)
save_name ='NR'+str(int(NR))+'_NS'+str(int(NS))+'_average'
ax.set_xlim(0, alpha0[-1]*(1.01))
ax.set_xlabel(r'$\alpha_0$')
    #ax.set_ylabel(r'$\max_{(\gamma_0, S_0)\in [0,1]^2}|\langle$Re($\lambda_1$)$\rangle|$')
ax.set_ylabel(r'$\max|\langle$Re($\lambda_1$)$\rangle|$')
ax.set_yscale('linear')
ax.ticklabel_format(axis="both", style="sci", scilimits=(-2,2))
ax.legend()
fig.tight_layout()
fig.savefig('additional_plots/largest_eigenvalue_varying_syntrophy_'+save_name+'.pdf')
plt.close()


fig=plt.figure()
ax=fig.add_subplot(111)
for i in range(N_alphamodes):
    numb_mat=[]
    for j in range(len(alpha0)):
        numb_mat.append(len(cont_matrices[i][j])/Nmatrices)
    ax.plot(alpha0, numb_mat, label=label[i])
save_name ='_cont_mat_NR'+str(int(NR))+'_NS'+str(int(NS))+'_average'
ax.set_xlim(0, alpha0[-1]*(1.01))
ax.set_xlabel(r'$\alpha_0$')
ax.set_ylabel(r'Prop. feasible matrices')
ax.set_yscale('linear')
ax.ticklabel_format(axis="both", style="sci", scilimits=(-2,2))
ax.legend()
fig.tight_layout()
fig.savefig('additional_plots/largest_eigenvalue_varying_syntrophy_'+save_name+'.pdf')
plt.close()

fig=plt.figure()
ax=fig.add_subplot(111)
for i in range(N_alphamodes):
    av_conn=[]
    for j in range(len(alpha0)):
        av_conn.append(np.mean(connectance[cont_matrices[i][j]]))
    ax.plot(alpha0, av_conn, label=label[i])
save_name ='_conn_cont_mat_NR'+str(int(NR))+'_NS'+str(int(NS))+'_average'
ax.set_xlim(0, alpha0[-1]*(1.01))
ax.set_xlabel(r'$\alpha_0$')
ax.set_ylabel(r'Av. conn. of feasible $G$')
ax.set_yscale('linear')
ax.ticklabel_format(axis="both", style="sci", scilimits=(-2,2))
ax.legend()
fig.tight_layout()
fig.savefig('additional_plots/largest_eigenvalue_varying_syntrophy_'+save_name+'.pdf')
plt.close()

fig=plt.figure()
ax=fig.add_subplot(111)
for i in range(N_alphamodes):
    av_nest=[]
    for j in range(len(alpha0)):
        av_nest.append(np.mean(nestedness[cont_matrices[i][j]]))
    ax.plot(alpha0, av_nest, label=label[i])
save_name ='_nest_cont_mat_NR'+str(int(NR))+'_NS'+str(int(NS))+'_average'
ax.set_xlim(0, alpha0[-1]*(1.01))
ax.set_xlabel(r'$\alpha_0$')
ax.set_ylabel(r'Av. nest. of feasible $G$')
ax.set_yscale('linear')
ax.ticklabel_format(axis="both", style="sci", scilimits=(-2,2))
ax.legend()
fig.tight_layout()
fig.savefig('additional_plots/largest_eigenvalue_varying_syntrophy_'+save_name+'.pdf')
plt.close()
