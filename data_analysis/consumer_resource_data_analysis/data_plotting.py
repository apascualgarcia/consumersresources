import numpy as np
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import math
import matplotlib.tri as tr

#======= GENERAL PLOTTING VARIABLES ========#
np.set_printoptions(threshold=sys.maxsize)

mpl.rcParams['lines.linewidth']=1.5
mpl.rcParams['lines.markersize']=12
mpl.rcParams['lines.markeredgewidth']=3
mpl.rcParams['lines.marker']='.'
mpl.rcParams['lines.linestyle']='solid'
mpl.rcParams['xtick.labelsize']=15
mpl.rcParams['ytick.labelsize']=15

alpha0=np.array([0., 0.0013, 0.0026, 0.0039, 0.0052, 0.0065, 0.0078, 0.0091, 0.0104, 0.014])
all_nestedness=[0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]
all_connectance=[0.08, 0.13, 0.18, 0.23, 0.28, 0.33, 0.38, 0.43]

alpha_mode=['fully_connected','random_structure', 'optimized_matrix', 'binary_metamatrix']
alpha_mode_colours=dict({'fully_connected':'blue', 'binary_metamatrix':'orange', 'optimized_matrix':'green', 'random_structure':'red'})
alpha_mode_sym = dict({'fully_connected': 'P', 'binary_metamatrix': 'D', 'optimized_matrix':'o', 'random_structure': 'd'})

alpha_mode_label=dict({'fully_connected':'Fully Connected', 'binary_metamatrix': 'Binary Optimized', 'optimized_matrix': 'Optimized', 'random_structure':'Random'})
legend_titles=dict({'nestG': 'Consumption \n overlap '+r'$\eta_\mathrm{G}$', 'connG': 'Consumption \n connectance '+r'$\kappa_\mathrm{G}$'})

labels=dict({'nestG': r'Consumption overlap $\eta_\mathrm{G}$', 'connG': r'Consumption connectance',
            'E': r'Objective function', 'connA': r'Syntrophy connectance',
            'nestA': r'Syntrophy overlap', 'feasible decay rate': r'Feasibility decay',
            'ld stable decay rate': r'Dynamical stability decay', 'alpha0': r'Syntrophy strength $\alpha_0$',
            'av. dominant eigenvalue':r'Mean rate of return', 'feasible volume':r'Feasible volume', 'dynamically stable volume':r'Dynamically stable volume'})

nest_colours=[plt.cm.get_cmap('jet_r')(i/len(all_nestedness)) for i in range(len(all_nestedness))]
conn_colours=[plt.cm.get_cmap('jet_r')(i/len(all_connectance)) for i in range(len(all_connectance))]

cmap = plt.cm.get_cmap('jet_r')

#======== DO NOT MODIFY BELOW ===============#
import pandas as pd
import consumer_resource_data_analysis.math_functions as mf
import consumer_resource_data_analysis.data_stats as stats
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch
import matplotlib as mpl


def plot_figure_2A(ax, df, plotting_dict):
    ax = plot_volumes(ax, df, plotting_dict['width'], [plotting_dict['intershift'],plotting_dict['intrashift']], plotting_dict['alpha0'], plotting_dict['alpha mode'],data_type="feasible")
    ax.set_ylabel(labels['feasible volume'])
    ax.set_yscale('linear')
    return ax


def plot_figure_2B(ax, df, plotting_dict):
    ax = plot_volumes(ax, df, plotting_dict['width'], [plotting_dict['intershift'],plotting_dict['intrashift']], plotting_dict['alpha0'], plotting_dict['alpha mode'],data_type="dynamically stable")
    ax.set_ylabel(labels['dynamically stable volume'])
    ax.set_yscale('linear')
    return ax


def plot_figure_2C(ax, df, plotting_dict):
    ax = plot_volumes(ax, df, plotting_dict['width'], [plotting_dict['intershift'],plotting_dict['intrashift']], plotting_dict['alpha0'], plotting_dict['alpha mode'],data_type="av. dominant eigenvalue")
    ax.set_ylabel(r'Mean rate of return')
    ax.set_yscale('log')
    return ax

def plot_figure_2E(ax, df, plotting_dict):
    ax = plot_p_values_vs_alpha0(ax, df, plotting_dict, 'av. dominant eigenvalue')
    return ax

def plot_figure_2F(ax, df, plotting_dict):
    ax = plot_p_values_vs_alpha0(ax, df, plotting_dict, 'feasible volume')
    return ax

def plot_figure_2G(ax, df, plotting_dict):
    ax = plot_p_values_vs_alpha0(ax, df, plotting_dict, 'dynamically stable volume')
    return ax

def plot_p_values_vs_alpha0(axes, data_frame, plotting_dict, data_type):
    # take all possibles alpha_modes pairs
    df = data_frame.copy()
    a0modes = plotting_dict['alpha mode']
    alpha0s = plotting_dict['alpha0']
    pairs = [(a0modes[i], a0modes[j]) for i in range(len(a0modes)) for j in range(i+1, len(a0modes))]
    for pair in pairs:
        a0mode1 = pair[0]
        a0mode2 = pair[1]
        p_vals = []
        a_vals = []
        for a0 in alpha0s:
            distrib1 = df[(df['alpha0']==a0) & (df['A-mode']==a0mode1)][data_type].to_numpy()
            distrib2 = df[(df['alpha0']==a0) & (df['A-mode']==a0mode2)][data_type].to_numpy()

            pval = stats.pvalue(distrib1, distrib2)

            if not(distrib1.any()) or not(distrib2.any()):
                print("One of the distributions is empty for alpha0=",a0," :")
                print(" A-mode = ", a0mode1," distribution = ", distrib1)
                print(" A-mode = ", a0mode2," distribution = ", distrib2)

            p_vals.append(pval)
            a_vals.append(a0)
        axes.plot(a_vals, p_vals, label=alpha_mode_label[a0mode1]+'-'+alpha_mode_label[a0mode2])
        axes.legend()
        axes.set_xlabel(labels['alpha0'])
        axes.set_ylabel(r'$p$-value (Mann-Whitney test)')
        axes.set_title(labels[data_type])

    return axes

def plot_volumes(ax, data_frame_, width, shift, alpha0_, alpha_mode_, data_type):

    df = data_frame_.copy()

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

def plot_against_matrix_properties(axs, data_frame_, plotting_properties, data_type):
    data = data_frame_.copy()
    ms = 8
    amodes = plotting_properties['alpha mode']
    for j in range(len(amodes)):
        amode = amodes[j]
        alpha0=0.0013
        for i in range(len(all_connectance)):
            connG = all_connectance[i]
            indices = [x for x in range(data.shape[0]) if mf.closest_element_in_list(data.loc[x]['G connectance'], all_connectance)==connG and data.loc[x]['A-mode']==amode and data.loc[x]['alpha0']==alpha0]
            to_plot = data.loc[indices]
            if data_type == 'av. dominant eigenvalue':
                to_plot[data_type]=-to_plot[data_type]
            to_plot = to_plot.sort_values(by='G nestedness')
            to_plot.plot(x = 'G nestedness', y= data_type, ax=axs[j][1], label=r'$'+str(connG)+'$',
                    linestyle='solid', color=conn_colours[i], marker=alpha_mode_sym[amode], markersize=ms)
        min_nest = np.min(data['G nestedness'])
        max_nest = np.max(data['G nestedness'])
        axs[j][1].set_xlim(min_nest-0.1*(max_nest-min_nest), max_nest+0.1*(max_nest-min_nest))
        axs[j][1].set_xlabel(labels['nestG'])
        axs[j][1].set_ylabel(labels[data_type])
        axs[j][1].legend(bbox_to_anchor=(1.,1.), loc='upper left', title=legend_titles['connG'], fontsize=13, title_fontsize=14)
        axs[j][1].set_title(alpha_mode_label[amode])

    for j in range(len(amodes)):
        amode = amodes[j]
        alpha0=0.0013
        for i in range(len(all_nestedness)):
            nestG = all_nestedness[i]
            indices = [x for x in range(data.shape[0]) if mf.closest_element_in_list(data.loc[x]['G nestedness'], all_nestedness)==nestG and data.loc[x]['A-mode']==amode and data.loc[x]['alpha0']==alpha0]
            to_plot = data.loc[indices]
            if data_type == 'av. dominant eigenvalue':
                to_plot[data_type]=-to_plot[data_type]
            to_plot = to_plot.sort_values(by='G connectance')
            to_plot.plot(x = 'G connectance', y= data_type, ax=axs[j][0], label=r'$'+str(nestG)+'$',
                    linestyle='solid', color=nest_colours[i], marker=alpha_mode_sym[amode], markersize=ms)
        min_conn = np.min(data['G connectance'])
        max_conn = np.max(data['G connectance'])
        axs[j][0].set_xlim(min_conn-0.1*(max_conn-min_conn), max_conn+0.1*(max_conn-min_conn))
        axs[j][0].set_xlabel(labels['connG'])
        axs[j][0].set_ylabel(labels[data_type])
        axs[j][0].legend(bbox_to_anchor=(1.,1.), loc='upper left', title=legend_titles['nestG'], fontsize=13, title_fontsize=14)
        axs[j][0].set_title(alpha_mode_label[amode])


    return axs

def plot_levels(data, colors_, level_name, alphamodes):
    # On récupère la liste des alpha0
    alpha0s = sorted(pd.unique(data['alpha0']))
    Nalpha0s = len(alpha0s)

    # Création des figures - canevas vide
    N_alphamodes = len(alphamodes)
    if N_alphamodes>1:
        fig, axs = plt.subplots(1, N_alphamodes, sharey=True, sharex=True, figsize=(3.5*N_alphamodes,4.5))
    else:
        fig = plt.figure()
        axs = fig.add_subplot(111)

    # find max alpha0 for which there is a non-zero fully feasible/stable/etc. volume -> this is the goal of the to_plot variable
    gamma0s = sorted(pd.unique(data['gamma0']))
    S0s = sorted(pd.unique(data['S0']))

    i=0
    for alphamode in alphamodes:
        # data analysis, find the right data subset
        to_analyse = data[data['alpha_mode']==alphamode][['gamma0', 'S0', 'alpha0', level_name]].sort_values(['gamma0', 'S0','alpha0'])
        max_ff = []
        for gamma0 in gamma0s:
            for S0 in S0s:
                sub_data = to_analyse[(to_analyse['gamma0']==gamma0)&(to_analyse['S0']==S0)]
                alpha0_max = np.max(sub_data[sub_data[level_name]==1.]['alpha0'])
                max_ff.append([gamma0, S0, alpha0_max])
        max_ff = np.array(max_ff)
        triang = tr.Triangulation(max_ff[:,0], max_ff[:,1])
        isbad = np.isnan(max_ff[:,2])
        # masque les valeurs nan
        mask = np.any(np.where(isbad[triang.triangles], True, False), axis=1)
        triang.set_mask(mask)

        axs[i].set_aspect('equal')
        im = axs[i].tricontourf(triang, max_ff[:,2], colors=colors_)
        axs[i].set_xlabel(r'$\gamma_0$')
        axs[i].set_title(alpha_mode_label[alphamode])
        axs[i].set_xticks([0.01, 0.5, 1])
        axs[i].set_xticklabels([0.01, 0.5, 1])

        i+=1


    axs[0].set_ylabel(r'$S_0$')
    axs[0].set_yticks([0.01, 0.5, 1])
    axs[0].set_yticklabels([0.01, 0.5, 1])

    cbar_ax = fig.add_axes([0.125, 0.1375, 0.775, 0.02])

    cbar=add_discrete_color_bar(cbar_ax, colors_, alpha0, [str(a) for a in alpha0])
    cbar.set_label(r'$\alpha_0$')


    fig.subplots_adjust(bottom=0.2, top=1)

    return fig, axs

def add_discrete_color_bar(ax_,colors_, levels, tick_labels, N_alphamodes=3):
    cmap_ = mpl.colors.ListedColormap(colors_)
    cmap_.set_over('0.25')
    cmap_.set_under('0.75')

    bounds = [i for i in range(len(levels)+1)]
    norm = mpl.colors.BoundaryNorm(bounds, cmap_.N)
    cbar = mpl.colorbar.ColorbarBase(ax_, cmap=cmap_,
                                norm=norm,
                                extend='neither',
                                spacing='uniform',
                                orientation='horizontal')


    cbar.set_ticks(ticks=[a +0.5 for a in bounds])
    cbar.set_ticklabels(tick_labels)

    return cbar

def plot_parameters_fixed_alpha0(data, level_name, alpha0_val, alphamodes, cmap = mpl.cm.get_cmap('bwr_r')):

    # Création des figures - canevas vide
    N_alphamodes = len(alphamodes)
    if N_alphamodes>1:
        fig, axs = plt.subplots(1, N_alphamodes, sharey=True, sharex=True, figsize=(3.5*N_alphamodes,4.5))
    else:
        fig = plt.figure()
        axs = [fig.add_subplot(111)]


    i=0
    for alphamode in alphamodes:
        to_plot = data[(data['alpha_mode']==alphamode) & (data['alpha0']==alpha0_val)][['gamma0', 'S0',level_name]].to_numpy()
        triang = tr.Triangulation(to_plot[:,0], to_plot[:,1])
        isbad = np.isnan(to_plot[:,2])
        # masque les valeurs nan
        mask = np.any(np.where(isbad[triang.triangles], True, False), axis=1)
        triang.set_mask(mask)

        axs[i].set_aspect('equal')
        im = axs[i].tricontourf(triang, to_plot[:,2], cmap=cmap)
        axs[i].set_xlabel(r'$\gamma_0$')
        #axs[i].set_title(alpha_mode_label[alphamode])
        axs[i].set_xticks([0.01, 0.5, 1])
        axs[i].set_xticklabels([0.01, 0.5, 1])

        i+=1


    axs[0].set_ylabel(r'$S_0$')
    axs[0].set_yticks([0.01, 0.5, 1])
    axs[0].set_yticklabels([0.01, 0.5, 1])

    norm= mpl.colors.Normalize(vmin=0, vmax=1)
    sm = plt.cm.ScalarMappable(norm=norm, cmap = cmap)
    sm.set_array([])

    cbar = fig.colorbar(sm, ticks=[0,0.25, 0.5, 0.75,1.])

    return fig, axs
