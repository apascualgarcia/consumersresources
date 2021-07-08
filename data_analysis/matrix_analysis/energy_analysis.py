import common_features.mpl_params
from common_features.functions import asymptote
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#### FUNCTION DEFINITION NOT CUSTOMIZABLE PART #########

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w
def plot_moving_average(ax, x_axis, y_axis, symbol, legend):
    #ax.plot(x_axis, y_axis, 'o', markersize=0.5)
    ax.plot(moving_average(y_axis, 15000), marker='', linestyle='solid', linewidth=2, label=legend)
    ax.set_ylabel(symbol)
    ax.set_xlabel(r'steps')
    return

# plots energy, nestedness and connectance on the axis that are provided
def plot_data_optimized(axis, min_data, label=''):
    energy = min_data[:,0]
    nest = min_data[:,1]
    conn = min_data[:, 2]
    T = min_data[:, 3]
    steps = np.linspace(start=0, stop=len(T)-1, num=len(T))

    plot_moving_average(axis[0], steps, energy, r'$\langle E \rangle$', legend=label)
    plot_moving_average(axis[1], steps, nest,r'$\langle \eta_A\rangle$', legend=label)
    plot_moving_average(axis[2], steps, conn, r'$\langle \kappa_A\rangle$', legend=label)

    return

###### END OF NON CUSTOMIZABLE PART : START OF THE SCRIPT #########


umatrix_data_path = "data_output/RandTrix_Nr25_Nc25_Nest0.45_Conn0.424_optimal_alpha.txt"
modes = [
    #"_alpha0=0.5_intra_specific_syntrophy=allowed_verbose-level=1_gamma0=1_",
    #"_alpha0=0.5_intra_specific_syntrophy=not_allowed_verbose-level=1_gamma0=1_",
    "_alpha0=1_intra_specific_syntrophy=allowed_verbose-level=1_gamma0=1_",
    "_alpha0=1_intra_specific_syntrophy=not_allowed_verbose-level=1_gamma0=1_"
    ]

modes_label = [
    #r"$\alpha_0 = 0.5$ ISS",
    #r"$\alpha_0 = 0.5$ no ISS",
    r"$\alpha_0 = 1$ ISS",
    r"$\alpha_0 = 1$ no ISS"
    ]

save_path='plots/RandTrix_Nr25_Nc25_Nest0.4_Conn0.1232_optimal_alpha'
g_nest='0.45'
g_conn='0.424'
title = r'$\kappa_G ='+g_conn+', \ \eta_G='+g_nest+'$'
figs_save_name = ["_energy.png", "_nestedness.png", "_connectance.png"]

figs = []
axis = []

for i in range(3):
    figs.append(plt.figure(i))
    axis.append(figs[i].add_subplot(111))
    axis[i].set_title(title)

for i in range(len(modes)):
    plot_data_optimized(axis, min_data=np.loadtxt(umatrix_data_path+modes[i]+'energy'), label=modes_label[i])

for i in range(3):
    axis[i].legend()
    figs[i].tight_layout()
    figs[i].savefig(save_path+figs_save_name[i], dpi=200)
