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
mat_index=4
mat_index_2=23
exp_coeff=[0.2,0.25,0.111111,0.111111,0.125,0.0833333,0.0833333,0.0833333,0.1,0.111111,0.0769231,0.0666667,0.0666667,0.0769231,0.0769231,0.0714286,0.0769231,0.0588235,0.0588235,0.0625,0.0666667,0.0625,0.0666667,0.0555556,0.0555556,0.0526316]

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

title_plot=r'$N_R='+str(int(NR[mat_index]))
title_plot+=r', N_S='+str(int(NS[mat_index]))+'$'


#### PLOT A COUPLE OF COOL MATRICES #####
#get critical points for fit
fit_indices=np.less(feasability[mat_index],0.6) & np.greater(feasability[mat_index], 0.4)
g0_to_fit=gamma0[mat_index][fit_indices]
S0_to_fit=S0[mat_index][fit_indices]

popt, pcov = curve_fit(func_fit, g0_to_fit, S0_to_fit)
print(popt, pcov)
fit_g0 = np.linspace(0.00, 1., 10000)
fit_S0 = np.ma.array([func_fit(x, *popt) for x in fit_g0])
fit_S0 = np.ma.masked_where(fit_S0 > 1., fit_S0)

#get critical points for fit
fit_indices_2=np.less(feasability[mat_index_2],0.6) & np.greater(feasability[mat_index_2], 0.4)
g0_to_fit_2=gamma0[mat_index_2][fit_indices_2]
S0_to_fit_2=S0[mat_index_2][fit_indices_2]

popt, pcov = curve_fit(func_fit, g0_to_fit_2, S0_to_fit_2)
print(popt, pcov)
fit_g0_2 = np.linspace(0.00, 1., 10000)
fit_S0_2 = np.ma.array([func_fit(x, *popt) for x in fit_g0_2])
fit_S0_2 = np.ma.masked_where(fit_S0_2 > 1., fit_S0_2)

expected_S0 = np.ma.array([exp_coeff[mat_index]/x for x in fit_g0])
expected_S0 = np.ma.masked_where(expected_S0 > 1., expected_S0)

expected_S0_2 = np.ma.array([exp_coeff[mat_index_2]/x for x in fit_g0_2])
expected_S0_2 = np.ma.masked_where(expected_S0_2 > 1., expected_S0_2)

cmap = cm.get_cmap('bwr_r')
fig1 =plt.figure(1)
ax1 = fig1.add_subplot(121, aspect='equal')
im = ax1.tricontourf(gamma0[mat_index], S0[mat_index], feasability[mat_index], cmap=cmap, levels=[0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.00001])
#ax1.plot(fit_g0, fit_S0, color='black', marker='None',linestyle='-.', linewidth=2)
ax1.plot(fit_g0, expected_S0, color='black', marker='None', linestyle='--', linewidth=2)
ax1.set_xlim(0., 1.)
ax1.set_ylim(0., 1.)
ax1.set_xlabel(r'$\gamma_0$')
ax1.set_ylabel(r'$S_0$')
title_axis=r'$\eta='+str(nestedness[mat_index])
title_axis+=r', \kappa='+str(connectance[mat_index])+'$'
ax1.set_title(title_axis)


ax12 = fig1.add_subplot(122, aspect='equal')
im = ax12.tricontourf(gamma0[mat_index], S0[mat_index_2], feasability[mat_index_2], cmap=cmap, levels=[0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.00001])
#ax12.plot(fit_g0_2, fit_S0_2, color='black', marker='None',linestyle='-.', linewidth=2)
ax12.plot(fit_g0_2, expected_S0_2, color='black', marker='None', linestyle='--', linewidth=2)

ax12.set_xlim(0., 1.)
ax12.set_ylim(0., 1.)
ax12.set_yticklabels([])
ax12.set_xlabel(r'$\gamma_0$')
title_axis=r'$\eta='+str(nestedness[mat_index_2])
title_axis+=r', \kappa='+str(connectance[mat_index_2])+'$'
ax12.set_title(title_axis)

#fig1.subplots_adjust(right=3)
cbar=fig1.colorbar(im, ax=[ax1, ax12], orientation='horizontal', shrink=0.8, pad=-0.4)
cbar.set_label(r'$\mathcal{F}(\gamma_0, S_0, G)$')
fig1.suptitle(title_plot)
fig1.tight_layout()


# plot the common feasability volume
# ax1.tricontour(triang, in_feas_vol, colors='k', levels=2)


### COMMON FEASIBILITY PART #####
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
#mask = np.all(np.where(do_not_plot[triang.triangles], True, False), axis=1)
#triang.set_mask(mask)

expected_min_S0 = np.ma.array([0.0526316/x for x in fit_g0])
expected_min_S0 = np.ma.masked_where(expected_min_S0 > 1., expected_min_S0)

min_S0_index = [i for i in range(len(feasability[mat_index])-1) if (not(do_not_plot[i]) and do_not_plot[i+1])]
popt, pcov = curve_fit(func_fit, gamma0[mat_index][min_S0_index], S0[mat_index][min_S0_index])
print(popt, pcov)

fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111, aspect='equal')
ax2.plot(fit_g0, expected_min_S0, color='black', marker='None',linestyle='--', linewidth=2, label=r'Theoretical minimum $S_0$')
ax2.tricontourf(triang, in_feas_vol,colors=['red', 'blue'], label='Common feasible region', levels=[0,0.5,1.001])
ax2.set_xlim(0., 1.)
ax2.set_ylim(0., 1.)
ax2.set_xlabel(r'$\gamma_0$')
ax2.set_ylabel(r'$S_0$')

fig2.tight_layout()

fig1.savefig(plot_folder+'/typical_feasibility_volume.pdf')
fig2.savefig(plot_folder+'/common_feasibility_volume.pdf')

deviations_from_theory = []

for i in range(len(exp_coeff)):
    fit_indices=np.less(feasability[i],0.6) & np.greater(feasability[i], 0.4)
    g0_to_fit=gamma0[i][fit_indices]
    S0_to_fit=S0[i][fit_indices]

    popt, pcov = curve_fit(func_fit, g0_to_fit, S0_to_fit)
    diff = (popt[0]-exp_coeff[i])/popt[0]
    deviations_from_theory.append(diff)
print([(connectance[i],deviations_from_theory[i]) for i in range(len(exp_coeff))])
rcParams.update({'font.size': 22})

fig3 = plt.figure(3)
ax3 = fig3.add_subplot(111)
nestedness_list=sorted(list(set(nestedness)))
for nest in nestedness_list:
    ind_to_plot=[i for i in range(len(exp_coeff)) if nestedness[i]==nest]
    dev = [deviations_from_theory[i] for i in ind_to_plot]
    conn = [connectance[i] for i in ind_to_plot]
    ax3.plot(conn, dev, linestyle='', label=r'$\eta='+str(nest)+'$', markersize=7, markeredgewidth=6)
ax3.set_xlabel(r'$\kappa$')
ax3.set_ylabel(r'$\Delta_G$')
ax3.legend(bbox_to_anchor=(1.0, 1.3), fontsize=17)

fig4 = plt.figure(4)
ax4 = fig4.add_subplot(111)
connectance_list=sorted(list(set(connectance)))
for conn in connectance_list:
    ind_to_plot=[i for i in range(len(exp_coeff)) if connectance[i]==conn]
    dev = [deviations_from_theory[i] for i in ind_to_plot]
    nest = [nestedness[i] for i in ind_to_plot]
    ax4.plot(nest, dev, linestyle='',label=r'$\kappa='+str(conn)+'$', markersize=7, markeredgewidth=6)
ax4.set_xlabel(r'$\eta$')
ax4.set_ylabel(r'$\Delta_G$')
ax4.legend(bbox_to_anchor=(1.0, 1.0))

fig3.tight_layout()
fig4.tight_layout()

fig3.savefig(plot_folder+'/feasibility_away_from_theory_fixed_nestedness.pdf')
fig4.savefig(plot_folder+'/feasibility_away_from_theory_fixed_connectance.pdf')

plt.show()
