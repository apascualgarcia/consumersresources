import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy.optimize import curve_fit

save_folder='plots/Phase_diagrams/Fully_connected_network'
data_folder = 'data_output'
filename = 'phase_diagram_configuration_comparison_NR25_NS25_s05_a05'
title = '$\alpha_0 = 0.5, \ \sigma_0 = 0.5$'


def power_fit(x, a, b):
    return b*np.power(x-1,a)


data = np.loadtxt(data_folder+'/'+filename+'.out')
NR = data[:, 0]
NS = data[:, 1]

# stability matrix has NR as rows and NS as columns
stability_matrix = []
for i in range(0, int(max(NR))):
    stability_matrix.append([])
    for j in range(0, int(max(NS))):
        index = data[:, 0:2].tolist().index([i+1,j+1])
        stability_matrix[i].append(data[index, 2:5])
stability_matrix = np.array(stability_matrix)

Nrange = np.linspace(1, min(len(stability_matrix), len(stability_matrix[0])), min(len(stability_matrix), len(stability_matrix[0])), dtype=int)
equi_unstable = [stability_matrix[i-1,i-1,0] for i in Nrange]
equi_marginal = [stability_matrix[i-1,i-1,1] for i in Nrange]
equi_stable = [stability_matrix[i-1,i-1,2] for i in Nrange]

# fitting to find critical exponents
# popt, trash = curve_fit(power_fit, Nrange, equi_unstable)
# opta,optb = popt

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
#ax1.plot(Nrange, [power_fit(x,opta, optb) for x in Nrange], linestyle='-', marker=None)
ax1.plot(Nrange, equi_unstable , label='Unstable')
ax1.plot(Nrange, equi_marginal, label='Marginal')
ax1.plot(Nrange, equi_stable, label='Stable')
ax1.legend()
ax1.set_xlabel(r'$N_R=N_S$')
ax1.set_ylabel(r'Proportion')
fig1.tight_layout()
fig1.savefig(save_folder+'/equispeciesresources_'+filename)

fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111, projection='3d')
ax2.plot_trisurf(NR, NS, data[:, 3])
ax2.view_init(azim=45, elev=26)
ax2.set_xlabel(r'$N_R$')
ax2.set_ylabel(r'$N_S$')
ax2.set_zlabel('Marginally stable systems')
fig2.tight_layout()
fig2.savefig(save_folder+'/surface_marginal_'+filename)

fig3 = plt.figure(3)
ax3 = fig3.add_subplot(111)
ax3.imshow(stability_matrix, origin='lower', interpolation='none', extent = (0.5, int(max(NS))+0.5, 0.5, int(max(NR))+0.5))
ax3.set_xlabel(r'$N_S$')
ax3.set_ylabel(r'$N_R$')
fig3.tight_layout()
fig3.savefig(save_folder+'/color_surf_'+filename)
