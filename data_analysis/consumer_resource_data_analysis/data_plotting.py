import numpy as np
import sys
import matplotlib as mpl

#======= GENERAL PLOTTING VARIABLES ========#
np.set_printoptions(threshold=sys.maxsize)

mpl.rcParams['lines.linewidth']=1.5
mpl.rcParams['lines.markersize']=12
mpl.rcParams['lines.markeredgewidth']=3
mpl.rcParams['lines.marker']='.'
mpl.rcParams['lines.linestyle']='solid'

alpha0=[0, 1.3e-3, 2.6e-3, 3.9e-3, 5.2e-3, 6.5e-3, 7.8e-3, 9.1e-3, 1.04e-2, 1.4e-2]
alpha_mode=['random_structure','optimized_matrix']
alpha_mode_colours=dict({'fully_connected':'blue', 'no_release_when_eat':'orange', 'optimized_matrix':'green', 'random_structure':'red'})
alpha_mode_sym = dict({'fully_connected': 'P', 'no_release_when_eat': 'D', 'optimized_matrix':'o', 'random_structure': 'd'})

alpha_mode_label=dict({'fully_connected':'Fully Connected', 'no_release_when_eat': 'NIS', 'optimized_matrix': 'Optimized', 'random_structure':'Random'})
legend_titles=dict({'nestG': 'Consumption \n overlap '+r'$\eta_\mathrm{G}$', 'connG': 'Consumption \n connectance '+r'$\kappa_\mathrm{G}$'})

labels=dict({'nestG': r'Consumption overlap $\eta_\mathrm{G}$', 'connG': r'Consumption connectance $\kappa_\mathrm{G}$',
            'E': r'Objective function $ E $', 'connA': r'Syntrophy connectance $\kappa_\mathrm{A}$',
            'nestA': r'Syntrophy overlap $\eta_\mathrm{A}$', 'feasible decay rate': r'Feasibility decay',
            'ld stable decay rate': r'Dynamical stability decay', 'alpha0': r'Syntrophy strength $\alpha_0$'})

#======== DO NOT MODIFY BELOW ===============#
import pandas as pd
import consumer_resource_data_analysis.math_functions as mf
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch
import matplotlib as mpl
import matplotlib.pyplot as plt


def plot_figure_2A(ax, df, plotting_dict):
    ax = plot_volumes(ax, df, plotting_dict['width'], [plotting_dict['intershift'],plotting_dict['intrashift']], plotting_dict['alpha0'], plotting_dict['alpha mode'],data_type="feasible")
    ax.set_ylabel(r'Feasible volume')
    ax.set_yscale('linear')
    return ax


def plot_figure_2B(ax, df, plotting_dict):
    ax = plot_volumes(ax, df, plotting_dict['width'], [plotting_dict['intershift'],plotting_dict['intrashift']], plotting_dict['alpha0'], plotting_dict['alpha mode'],data_type="dynamically stable")
    ax.set_ylabel(r'Dynamically stable volume')
    ax.set_yscale('linear')
    return ax


def plot_figure_2C(ax, df, plotting_dict):
    ax = plot_volumes(ax, df, plotting_dict['width'], [plotting_dict['intershift'],plotting_dict['intrashift']], plotting_dict['alpha0'], plotting_dict['alpha mode'],data_type="av. dominant eigenvalue")
    ax.set_ylabel(r'Mean rate of return')
    ax.set_yscale('log')
    return ax

def plot_volumes(ax, df, width, shift, alpha0_, alpha_mode_, data_type):
    intershift = shift[0]
    intrashift = shift[1]

    N_alphamodes = len(alpha_mode_)
    N_alpha0 = len(alpha0_)
    L = N_alphamodes*(width+intrashift)-intrashift+intershift
    xs = 0.5*(width+intershift)

    for a in df['alpha0']:
        df['alpha0'] = df['alpha0'].replace([a], mf.closest_element_in_list(a, alpha0_))

    ticks=[]
    for i in range(N_alpha0):
        box_start = i*L
        box_end = box_start+L
        if i%2==0:
            ax.axvspan(xmin=box_start, xmax=box_end, color='black', alpha=0.1)
        xloc = 0.5*(box_end-box_start)+box_start
        ticks.append(xloc)

    legend_els=[]
    all_means = []
    all_positions = []
    for i in range(N_alphamodes):
        amode=alpha_mode_[i]
        data=df[df['A-mode']==amode]
        col = alpha_mode_colours[amode]
        marker = alpha_mode_sym[amode]
        facecol = to_rgb(col)+(0.5,)
        boxprops = dict(edgecolor=facecol, facecolor=facecol, linewidth=0)
        meanprops=dict(color=col, marker=marker, markeredgecolor=col, markerfacecolor=col, linestyle='solid', markersize=10)
        medianprops=dict(color=col, marker='')
        whiskerprops=dict(color=facecol, marker='')
        capprops=dict(color=facecol, marker='')
        flierprops=dict(marker=marker, markeredgecolor=col, markerfacecolor=col, markeredgewidth=1, markersize=5)
        means=[]
        to_plot=[]
        for a0 in alpha0_:
            string = data_type+' volume'
            if data_type == "av. dominant eigenvalue":
                string = 'av. dominant eigenvalue'
            volumes = data[data['alpha0']==a0][string].to_numpy()
            if data_type == "av. dominant eigenvalue":
                volumes=-volumes
            volumes = volumes[~np.isnan(volumes)]
            means.append(np.mean(volumes))
            to_plot.append(list(volumes))
        all_means.append(means)
        positions = np.linspace(start=xs+i*(width+intrashift), stop=xs+(N_alpha0-1)*L+i*(width+intrashift), num=N_alpha0)
        all_positions.append(positions)
        ax.boxplot(to_plot, positions=positions, widths=width,
                patch_artist=True , boxprops=boxprops,
                meanprops=meanprops, medianprops=medianprops,
                flierprops=flierprops, whiskerprops=whiskerprops,
                capprops=capprops)
        legend_els.append(Patch(facecolor=facecol, edgecolor=facecol, linewidth=0, label=alpha_mode_label[amode]))

    # Remove comment if needed to plot solid lines between mean values
    for i in range(N_alphamodes):
        amode = alpha_mode_[i]
        col = alpha_mode_colours[amode]
        marker = alpha_mode_sym[amode]
        ax.plot(all_positions[i], all_means[i], marker=marker, markerfacecolor='black', markersize=5,
                linestyle='None', color=col, markeredgecolor='black', markeredgewidth=0.5)

    ax.set_xticks(ticks=ticks)
    ax.set_xticklabels([a for a in alpha0_], rotation=45)
    ax.set_xlim(0, xs+(N_alpha0-1)*L+(N_alphamodes-1)*(width+intrashift)+(width+intershift)*0.5)
    ax.legend(handles=legend_els, bbox_to_anchor=(0., 1.02, 1., .102), loc='center',
           ncol=len(legend_els), borderaxespad=0., fontsize=12)
    ax.set_title('')
    ax.set_xlabel(labels['alpha0'])


    return ax
