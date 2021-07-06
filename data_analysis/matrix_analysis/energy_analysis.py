import common_features.mpl_params
from common_features.functions import asymptote
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w
def plot_moving_average(ax, x_axis, y_axis, symbol, legend):
    #ax.plot(x_axis, y_axis, 'o', markersize=0.5)
    ax.plot(moving_average(y_axis, 15000), marker='', linestyle='solid', linewidth=2, label=legend)
    ax.set_ylabel(symbol)
    ax.set_xlabel(r'steps')
    return


ccmatrix_data_path = "data_output/RandTrix_Nr25_Nc25_Nest0.45_Conn0.1856_optimal_alpha.txt_constant_connectance_energy"
umatrix_data_path = "data_output/RandTrix_Nr25_Nc25_Nest0.45_Conn0.1856_optimal_alpha.txt_unconstrained_energy"
save_path='plots/RandTrix_Nr25_Nc25_Nest0.45_Conn0.1856_optimal_alpha'

g_nest='0.45'
g_conn='0.19'


title = r'$\kappa_G ='+g_conn+', \ \eta_G='+g_nest+'$'

ccdata = np.loadtxt(ccmatrix_data_path)
ccenergy = ccdata[:,0]
ccnest = ccdata[:,1]
ccconn = ccdata[:, 2]
ccT = ccdata[:, 3]
ccsteps = np.linspace(start=0, stop=len(ccT)-1, num=len(ccT))

udata = np.loadtxt(umatrix_data_path)
uenergy = udata[:,0]
unest = udata[:,1]
uconn = udata[:, 2]
uT = udata[:, 3]
usteps = np.linspace(start=0, stop=len(uT)-1, num=len(uT))

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
ax1.set_title(title)
plot_moving_average(ax1, ccsteps, ccenergy, r'$\langle E \rangle$', r'constant $\kappa_A$')
plot_moving_average(ax1, usteps, uenergy, r'$\langle E \rangle$', r'unconstrained')
ax1.legend()
fig1.tight_layout()
fig1.savefig(save_path+"_energy.png", dpi=200)

print("Final energy, cc: ", asymptote(ccenergy, 2000)," unc: ", asymptote(uenergy, 2000) )



fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111)
ax2.set_title(title)
plot_moving_average(ax2, ccsteps, ccnest,r'$\langle \eta_A \rangle$', r'constant $\kappa_A$')
plot_moving_average(ax2, usteps, unest,r'$\langle \eta_A\rangle$', r'unconstrained')
ax2.legend()
fig2.tight_layout()
fig2.savefig(save_path+"_nestedness.png", dpi=200)
print("Final nestedness, cc: ", asymptote(ccnest, 2000)," unc: ", asymptote(unest, 2000) )




fig3 = plt.figure(3)
ax3 = fig3.add_subplot(111)
ax3.set_title(title)
plot_moving_average(ax3, ccsteps, ccconn, r'$\langle \kappa_A\rangle$', r'constant $\kappa_A$')
plot_moving_average(ax3, usteps, uconn, r'$\langle \kappa_A\rangle$', r'unconstrained')
ax3.legend()
fig3.tight_layout()
fig3.savefig(save_path+"_connectance.png", dpi=200)
print("Final connectance, cc: ", asymptote(ccconn, 2000)," unc: ", asymptote(uconn, 2000) )




#x4 = fig1.add_subplot(224)
#ax4.plot(T, 'o', markersize=0.5)
#ax4.set_yscale('log')
#ax4.set_xlabel(r'steps')
#ax4.set_ylabel(r'$T$(steps)')
