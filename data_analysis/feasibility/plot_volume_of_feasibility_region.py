import numpy as np
import matplotlib.pyplot as plt
from common_functions import remove_strings_from_file
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import sys
import copy
import matplotlib.tri as tr
from scipy.optimize import curve_fit
from matplotlib import rcParams

alpha_mode=['random_structure', 'no_release_when_eat', 'optimal_matrix']
label=['fully connected', 'no release when eat', 'optimal LRI']
filename = 'feasibility/common_feasibility_volume_NR25_NS25_full_rank_opt_consumption_mat_NR25_NS25'
alpha0=[0, 0.0013, 0.0026, 0.0039, 0.0052, 0.0065, 0.0078, 0.0091, 0.0104, 0.014]

cfv_region = []
fig=plt.figure()
ax = fig.add_subplot(111)
k=0
for al_mo in alpha_mode:
    local_vector=[]
    for a in alpha0:
        file = filename+'_'+al_mo+'_optimal_LRI_alpha0='+str(a)+'.out'
        local_data=np.loadtxt(file)
        local_vector.append(len(local_data))
    no_synt_vol = local_vector[0]
    for i in range(len(local_vector)):
        local_vector[i]=local_vector[i]/no_synt_vol
    cfv_region.append(local_vector)
    ax.plot(alpha0, local_vector, label=label[k])
    k+=1
ax.set_xlabel(r'$\alpha_0$')
ax.set_ylabel(r'Vol $\mathcal{V}^1(\alpha_0)$/Vol $\mathcal{V}^1(0)$')
ax.set_yscale('linear')
ax.legend()
fig.tight_layout()
fig.savefig('plots/measure_common_feasibility_volume_varying_syntrophy.pdf')

# fig2 = plt.figure()
# ax = fig2.add_subplot(111)
# ax.plot(alpha0, [cfv_region[0][i]-cfv_region[1][i] for i in range(len(cfv_region[0]))])
# ax.set_xlabel(r'$\alpha_0$')
# ax.set_ylabel(r'Vol $\mathcal{V}^1(\alpha_0)$/Vol $\mathcal{V}^1(0)$')
# fig2.tight_layout()
# plt.show()
