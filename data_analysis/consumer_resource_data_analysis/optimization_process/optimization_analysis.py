import consumer_resource_data_analysis as cf
from consumer_resource_data_analysis.data_plotting import labels
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#### FUNCTION DEFINITION NOT CUSTOMIZABLE PART #########
def plot_moving_average(ax, x_axis, y_axis, symbol, legend):
    average = 1
    ax.plot(x_axis, y_axis, 'o', markersize=0.5, label=legend)
    ax.set_xscale('linear')
    #ax.plot(x_axis, cf.moving_average(y_axis, average), marker='', linestyle='solid', linewidth=2, label=legend)
    ax.set_ylabel(symbol)
    return

# plots energy, nestedness and connectance on the axis that are provided
def plot_data_optimized(axis, data, max_index=5000, label=''):
    plot_moving_average(axis[0], data['step'][:max_index], data['E'][:max_index], labels['E'], legend=label)
    plot_moving_average(axis[1], data['step'][:max_index], data['nestA'][:max_index],labels['nestA'], legend=label)
    plot_moving_average(axis[2], data['step'][:max_index], data['connA'][:max_index], labels['connA'], legend=label)

    return

###### END OF NON CUSTOMIZABLE PART : START OF THE SCRIPT #########


results_path = "/Users/Shared/Master/Master_Thesis/results/optimize_matrices_core_1.txt"
columns = ['step', 'T', 'E','nestA','connA','nestG','connG','eff_comp']
save_path = "/Users/Shared/Master/Master_Thesis/main_figures/monte_carlo/optimization_process"

runs = ['0','1','2','3']

g_nest='0.55'
g_conn='0.23'
title = r'$\kappa_G ='+g_conn+', \ \eta_G='+g_nest+'$'
figs_save_name = ["_energy.png", "_nestedness.png", "_connectance.png"]

figs = []
axis = []

for i in range(3):
    figs.append(plt.figure(i))
    axis.append(figs[i].add_subplot(111))
    axis[i].set_title(title)

for i in range(len(runs)):
    filename = results_path[:-5:]+str(i)+".txt"
    plot_data_optimized(axis, data=pd.read_csv(filename, dtype = float, sep=";", names=columns), label='Run '+str(i+1))

for i in range(3):
    axis[i].set_xlabel(r'Number of steps')
    axis[i].legend()
    figs[i].tight_layout()
    figs[i].savefig(save_path+figs_save_name[i], dpi=200)


plt.show()
