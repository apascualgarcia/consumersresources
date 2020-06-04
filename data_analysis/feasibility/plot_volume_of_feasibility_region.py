import numpy as np
import matplotlib.pyplot as plt
from common_functions import remove_strings_from_file
import consumer_resource_data_analysis as cf
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import sys
import copy
import matplotlib.tri as tr
from scipy.optimize import curve_fit
from matplotlib import rcParams
from consumer_resource_data_analysis import alpha_mode, label, alpha_mode_colours, alpha0

filename = 'feasibility/all_mat_feasibility_NR25_NS25_100_points_full_rank_opt_consumption_mat_NR25_NS25'
optimal_LRI_folder='optimal_LRI_Nr25_Nc25'
consumption_matrix_folder='optimal_matrices/consumption/Nr25_Nc25'
matrix_set='S_{25}'
ylabel=r'Vol$\left(\mathcal{F}_1^{S_{25}}(\alpha_0)\right)$'
feasibility_region = cf.load_data_region(alpha_mode, alpha0, filename, optimal_LRI_folder)
alpha0=np.array(alpha0)


NR=int(feasibility_region[0,0,0,0])
NS=int(feasibility_region[0,0,0,1])

fig=plt.figure()
ax=fig.add_subplot(111)
for i in range(len(alpha_mode)):
    volume=cf.common_full_volume_as_function_of_alpha0(feasibility_region, i, alpha0)
    fit_func=cf.exponential_function
    take_zero_points=True
    fitted_alpha0, fitted_vol, estimated_alpha_crit, estimated_vol_zero_syntrophy, (estimated_decay_rate, err)=cf.fit_shrinkage_curve(alpha0, volume, fit_func, take_zero_points)
    print(alpha_mode[i], ', estimated_decay_rate: ', estimated_decay_rate, '+/-', err)
    #ax.plot(fitted_alpha0, fitted_vol, color=alpha_mode_colours[i], linestyle='solid', marker='', alpha=0.7)
    ax.plot(alpha0, volume, label=label[i], color=alpha_mode_colours[i])
ax.set_xlabel(r'$\alpha_0$')
ax.set_yscale('log')
ax.set_ylabel(ylabel)
ax.legend()
fig.tight_layout()
fig.savefig('plots/common_feasibility_volume_NR'+str(NR)+'_NS'+str(NS)+'.pdf')

# cfv_region = []
# fig=plt.figure()
# ax = fig.add_subplot(111)
# k=0
# for al_mo in alpha_mode:
#     local_vector=[]
#     for a in alpha0:
#         file = filename+'_'+al_mo+'_optimal_LRI_alpha0='+str(a)+'.out'
#         local_data=np.loadtxt(file)
#         local_vector.append(len(local_data)/900)
#     no_synt_vol = local_vector[0]
#     for i in range(len(local_vector)):
#         local_vector[i]=local_vector[i]/no_synt_vol
#     cfv_region.append(local_vector)
#     fitted_vol, popt, pcov = cf.fit_data(cf.exponential_function, alpha0, local_vector)
#     print(popt, pcov)
#     ax.plot(alpha0, local_vector, label=label[k],markersize=10, linewidth=2.5, markeredgewidth=3)
#     k+=1
# ax.set_xlabel(r'$\alpha_0$')
# ax.set_ylabel(r'Vol $\left(\mathcal{F}^{S_M}_1(\alpha_0)\right)$ (normalized)')
# ax.set_yscale('log')
# ax.legend()
# fig.tight_layout()
# fig.savefig('plots/measure_common_feasibility_volume_varying_syntrophy.pdf')

# fig2 = plt.figure()
# ax = fig2.add_subplot(111)
# ax.plot(alpha0, [cfv_region[0][i]-cfv_region[1][i] for i in range(len(cfv_region[0]))])
# ax.set_xlabel(r'$\alpha_0$')
# ax.set_ylabel(r'Vol $\mathcal{V}^1(\alpha_0)$/Vol $\mathcal{V}^1(0)$')
# fig2.tight_layout()
# plt.show()
