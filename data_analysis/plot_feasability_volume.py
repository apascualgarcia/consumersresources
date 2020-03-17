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


np.set_printoptions(threshold=sys.maxsize)

data_folder='./data_output'
plot_folder='./plots'
matrices_folder='optimal_matrices/consumption/Nr25_Nc25'
filename='feasability_NR25_NS25'
mat_index=0

filename = data_folder+'/'+filename
remove_strings_from_file(matrices_folder, filename)
data=np.loadtxt(filename+'_filtered.out')

def func_fit(x, a,b):
    return a/x+b


NR=data[:,0]
NS=data[:,1]
nestedness=data[:,2]
connectance=data[:,3]
gamma0=data[:, 4::3]
S0=data[:, 5::3]
feasability=data[:, 6::3]

# plot typical matrix feasability volume
title_axis=r'$N_R='+str(int(NR[mat_index]))
title_axis+=r', N_S='+str(int(NS[mat_index]))
title_axis+=r', \eta='+str(nestedness[mat_index])
title_axis+=r', \kappa='+str(connectance[mat_index])
title_axis+='$'

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111, aspect='equal')
cmap = cm.get_cmap('bwr_r')
im1 = ax1.tricontourf(gamma0[mat_index], S0[mat_index], feasability[mat_index], cmap=cmap, levels=[0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.00001])
cbar=fig1.colorbar(im1)
ax1.set_title(title_axis, pad=12)
ax1.set_xlim(0., 1.)
ax1.set_xlabel(r'$\gamma_0$')
ax1.set_ylabel(r'$S_0$')
cbar.set_label(r'Feasibility')


# find the common feasability range
full_feas_indices=[j for j in range(len(feasability[0])) if feasability[0,j]==1.]
for i in range(1,len(data)):
    oldffi=copy.deepcopy(full_feas_indices)
    full_feas_indices=[j for j in range(len(feasability[i])) if (feasability[i,j]==1. and j in oldffi)]
in_feas_vol=[]
for i in range(len(feasability[mat_index])):
    if(i in full_feas_indices):
        in_feas_vol.append(1.)
    else:
        in_feas_vol.append(0.)

do_not_plot = np.less_equal(in_feas_vol,0.)
triang = tr.Triangulation(gamma0[mat_index], S0[mat_index])
mask = np.all(np.where(do_not_plot[triang.triangles], True, False), axis=1)
triang.set_mask(mask)

#get critical points for fit
fit_indices=np.less(feasability[mat_index],0.6) & np.greater(feasability[mat_index], 0.4)
g0_to_fit=gamma0[mat_index][fit_indices]
S0_to_fit=S0[mat_index][fit_indices]

popt, pcov = curve_fit(func_fit, g0_to_fit, S0_to_fit)
print(popt, pcov)
fit_g0 = np.linspace(0.01, 1., 10000)
fit_S0 = np.ma.array([func_fit(x, *popt) for x in fit_g0])
fit_S0 = np.ma.masked_where(fit_S0 > 1., fit_S0)

expected_min_S0 = np.ma.array([0.0526316/x for x in fit_g0])
expected_min_S0 = np.ma.masked_where(expected_min_S0 > 1., expected_min_S0)

min_S0_index = [i for i in range(len(feasability[mat_index])-1) if (not(do_not_plot[i]) and do_not_plot[i+1])]
popt, pcov = curve_fit(func_fit, gamma0[mat_index][min_S0_index], S0[mat_index][min_S0_index])
print(popt, pcov)

# plot the common feasability volume
# ax1.tricontour(triang, in_feas_vol, colors='k', levels=2)
ax1.plot(fit_g0, expected_min_S0, color='black', marker='None',linestyle='-', linewidth=1, label=r'Theoretical minimum $S_0$')
ax1.plot(fit_g0, fit_S0, linestyle='-.', marker='None', color='black', linewidth=1, label=r'Fit $\sim 1/\gamma_0$')
ax1.tricontourf(triang, in_feas_vol, hatches=['xx'],cmap='bwr_r', alpha=0., label='Common feasible region')
#ax1.legend()

fig1.tight_layout()

fig1.savefig(plot_folder+'/typical_feasibility_volume.pdf')
plt.show()
