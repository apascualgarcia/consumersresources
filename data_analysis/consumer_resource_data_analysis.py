import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import wilcoxon
import matplotlib.tri as tr
import matplotlib as mpl
import copy
import sys
import re
import pandas as pd
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch

np.set_printoptions(threshold=sys.maxsize)
mpl.rcParams['lines.linewidth']=1.5
mpl.rcParams['lines.markersize']=12
mpl.rcParams['lines.markeredgewidth']=3
mpl.rcParams['lines.marker']='.'
mpl.rcParams['lines.linestyle']='solid'

alpha_mode=['fully_connected', 'no_release_when_eat', 'optimal_matrix', 'random_structure']
alpha_mode_colours=dict({'fully_connected':'blue', 'no_release_when_eat':'orange', 'optimal_matrix':'red', 'random_structure':'green'})
alpha_mode_sym = dict({'fully_connected': 'P', 'no_release_when_eat': 'D', 'optimal_matrix':'o', 'random_structure': 'd'})
alpha_mode_label=dict({'fully_connected':'FC', 'no_release_when_eat': 'NIS', 'optimal_matrix': 'OM', 'random_structure':'RS'})
alpha0=[0, 1.3e-3, 2.6e-3, 3.9e-3, 5.2e-3, 6.5e-3, 7.8e-3, 9.1e-3, 1.04e-2, 1.4e-2]
all_nestedness=[0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]
all_connectance=[0.08, 0.13, 0.18, 0.23, 0.28, 0.33, 0.38, 0.43]
min_gamma0, max_gamma0=0.01, 1
min_S0, max_S0=0.01, 1
labels=dict({'nestG': r'Consumption overlap $\eta_G$', 'connG': r'Consumption connectance $\kappa_G$',
            'E': r'Ecosystem energy $ E $', 'connA': r'Syntrophy connectance $\kappa_A$',
            'nestA': r'Syntrophy overlap $\eta_A$', 'feasible decay rate': r'Feasibility decay',
            'ld stable decay rate': r'Dynamical stability decay'})
nest_colours=[plt.cm.get_cmap('jet_r')(i/len(all_nestedness)) for i in range(len(all_nestedness))]
conn_colours=[plt.cm.get_cmap('jet_r')(i/len(all_connectance)) for i in range(len(all_connectance))]
N_alphamodes=len(alpha_mode)

# alpha_mode=['optimal_matrix']
# label=['Modified LRI']
# N_alphamodes=1

# careful, only returns the real proba, ie the possible
def proba_dyna_if_feas(f_region_, d_region_, alpha_mode_index_, alpha0_index_, mat_index_):
    feasible_=f_region_[alpha_mode_index_, alpha0_index_, mat_index_, 6::3]
    dyn_stab_=d_region_[alpha_mode_index_, alpha0_index_, mat_index_, 6::3]

    equal_points=0
    feasible_points=0

    arr_to_return_ = []
    for i in range(len(feasible_)):
        if feasible_[i] > 0:
            feasible_points+=1
            arr_to_return_.append(dyn_stab_[i]/feasible_[i])
            if(dyn_stab_[i]==feasible_[i]):
                equal_points+=1
    # if feasible_points>0:
    #     print('A percentage ', equal_points/feasible_points, ' are 100 percent dynamically stable if feasible')

    return np.array(arr_to_return_)

# computes the volume (weigthed) of dynamical if feasible
def vol_dyna_if_feas(f_region_, d_region, alpha_mode_index_, alpha0_index_, mat_index_):
    to_return_ = 0
    weight=proba_dyna_if_feas(f_region_, d_region, alpha_mode_index_, alpha0_index_, mat_index_)
    length_=len(weight)
    if length_==0:
        return np.nan
    else:
        return np.mean(weight)

# gives you the average gamma0 (among all matrices), for the wanted alpha mode and alpha index
def average_coordinate_all_matrices(region_, alpha_mode_index, alpha0_index, coordinate_index):
    region_wanted=region_[alpha_mode_index, alpha0_index]
    Nmatrices=len(region_wanted)
    av_gamma0_=[]

    for i in range(Nmatrices):
        to_add = average_coordinate(region_wanted[i], coordinate_index)
        if(to_add > 0):
            av_gamma0_.append(to_add)
    print('Only', len(av_gamma0_), 'matrices had a non-zero volume')
    print(av_gamma0_)
    return np.mean(av_gamma0_)

# gives you the average coordinate (0 for gamma0, 1 for S0) as defined in the script, data is the data for one matrix
# at a given alpha mode and alpha0
def average_coordinate(data, i):
    quantity=data[6::3]
    gamma0=data[4+i::3]
    to_return_=0.
    Npoints=len(quantity)
    to_divide=np.sum(quantity)
    is_there_no_one=True
    for i in range(Npoints):
        if(quantity[i]>0):
            to_return_+=quantity[i]*gamma0[i]
        if(quantity[i]==1):
            is_there_no_one=False
    if(to_divide>0):
        to_return_/=to_divide
    # if(is_there_no_one):
    #     print("There is no fully dynamically stable point for that regime")
    return to_return_

def average_gamma0_all_matrices(region_, alpha_mode_index, alpha0_index):
    return average_coordinate_all_matrices(region_, alpha_mode_index, alpha0_index, 0)

def average_S0_all_matrices(region_, alpha_mode_index, alpha0_index):
    return average_coordinate_all_matrices(region_, alpha_mode_index, alpha0_index, 1)


# tells the common full volume for a given alphamode and alpha0
# it is computed by adding every point in the region pondered by its smallest
# quantity among all the matrices
def common_full_volume(region_, alpha_mode_index, alpha0_index):
    quantity_=region_[alpha_mode_index, alpha0_index, :, 6::3]
    gamma0_=region_[alpha_mode_index,alpha0_index,:,4::3]
    S0_=region_[alpha_mode_index,alpha0_index,:,5::3]

    Nmatrices_=len(region_[alpha_mode_index,alpha0_index])
    Npoints_=len(quantity_[0])

    volume_=0.
    for i_ in range(Npoints_):
        volume_ += min(quantity_[:, i_])

    return volume_/Npoints_

def common_full_volume_as_function_of_alpha0(region_, alpha_mode_index, alpha0_):
    to_plot=[]
    for k_ in range(len(alpha0_)):
        to_plot+=[common_full_volume(region_, alpha_mode_index, k_)]
    return to_plot

def plot_common_full_volume_shrinkage(ax_,region_, alpha_mode_index, alpha0_):
    to_plot=[]
    for k_ in range(len(alpha0_)):
        to_plot+=[common_full_volume(region_, alpha_mode_index, k_)]
    ax_.plot(alpha0_, to_plot)
    return


# region contains all data for each alphamode, alpha0 and network
def plot_common_region(region, alpha_mode_, colors_, labels_):

    fig, axs = plt.subplots(1, N_alphamodes, sharey=True, sharex=True, figsize=(3.5*N_alphamodes,4.5))
    quantity=region[:,:,:,6::3]
    gamma0=region[:,:,:,4::3]
    S0=region[:,:,:,5::3]
    Nmatrices=len(region[0,0])

    # do that computation for each alpha_mode

    for k in range(N_alphamodes):
        if N_alphamodes > 1:
            ax=axs[k]
        else:
            ax=axs
        # contains the indices which have quantity 1 for all matrices at different alpha0
        f_indices=[]
        # do that for all alpha0
        for l in range(len(region[k])):
            Npoints=len(quantity[k,l,0])
            # find the indices that have quantity=1 for this alpha0 and this alpha mode
            full_indices=full_indices_in_data_set(quantity[k,l,0])
            for m in range(1,Nmatrices):
                oldfi=copy.deepcopy(full_indices)
                full_indices_current_mat=full_indices_in_data_set(quantity[k,l,m])
                full_indices=[j for j in range(Npoints) if j in oldfi and j in full_indices_current_mat]
            f_indices.append(full_indices)
            #print('Full indices with method', len(full_indices))
        data_levels=[]
        for i in range(Npoints):
            level=-1
            j=0
            exit=False
            while(not(exit) and j < len(region[k])):
                if i in f_indices[j]:
                    level+=1
                else:
                    exit=True
                j+=1
            data_levels.append(level)
        max_level=max(data_levels)
        levels=[i for i in range(max_level+2)]
        triang = tr.Triangulation(gamma0[k,0,0], S0[k,0,0])
        ax.set_aspect('equal')
        ax.set_xlabel(r'$\gamma_0$')
        im=ax.tricontourf(triang, data_levels, levels=levels, colors=colors_)
        ax.set_xlim(min_gamma0, max_gamma0)
        ax.set_ylim(min_S0, max_S0)
        ax.set_title(labels_[k])
    if N_alphamodes > 1:
        axs[0].set_ylabel(r'$S_0$')
        axs[0].set_yticks([0, 0.5, 1])
        axs[0].set_yticklabels([0, 0.5, 1])

        fig.subplots_adjust(bottom=0.2, top=0.95)
    else:
        axs.set_ylabel(r'$S_0$')
        axs.set_yticks([0, 0.5, 1])
        axs.set_yticklabels([0, 0.5, 1])

    return fig, axs, im, levels

def filter_data(alpha_mode_, alpha0_, filename_, optimal_LRI_folder_, consumption_matrix_folder_):
    for al_mo in alpha_mode_:
        for a in alpha0_:
            file=filename_+'_'+al_mo+'_'+optimal_LRI_folder_+'_alpha0='+str(a)
            remove_strings_from_file(consumption_matrix_folder_, file)
    return
def load_data_region(alpha_mode_, alpha0_, filename_, optimal_LRI_folder_, type_=np.float64):
    region=[]
    for al_mo in alpha_mode_:
        local_vector=[]
        for a in alpha0_:
            file=filename_+'_'+al_mo+'_'+optimal_LRI_folder_+'_alpha0='+str(a)+'_filtered.out'
            local_data=np.loadtxt(file, dtype=type_)
            print('Loading file', file, 'which contains the data of', len(local_data), 'matrices')
            local_vector.append(local_data)
        region.append(local_vector)
    region=np.array(region)
    print('Dimensions of tableau : ', len(region), 'x', len(region[0]),'x',len(region[0,0]),'x', len(region[0,0,0]))
    return region


# returns the local dynamical stability levels (first index is alpha_mode, second is point)
# data has lds[alpha_mode][alpha0]Çpoint
def levels_different_alpha0(data):
    lds_levels=[]
    for i in range(len(data)):
        lds_levels.append(data_levels(data[i]))
    return np.array(lds_levels)

# returns the indices of the elements equal to one in the data set
def full_indices_in_data_set(data_):
    return [j for j in range(len(data_)) if data_[j]==1.]

# plot levels, the level -1 is not plotted, it has two indices, first is alphamode, second is alpha0
# colors is the colors chosen for all alpha0
def plot_levels(data, colors, labels_):
    gamma0=data[0,0,4::3]
    S0=data[0,0,5::3]
    if N_alphamodes>1:
        fig, axs = plt.subplots(1, N_alphamodes, sharey=True, sharex=True, figsize=(3.5*N_alphamodes,4.5))
    else:
        fig = plt.figure()
        axs = fig.add_subplot(111)
    function_levels=levels_different_alpha0(data)
    max_level=np.amax(function_levels)
    levels=[i for i in range(max_level+2)]
    triang = tr.Triangulation(gamma0, S0)
    for i in range(len(data)):
        to_plot = function_levels[i]
        if N_alphamodes > 1:
            axs[i].set_aspect('equal')
            im = axs[i].tricontourf(triang, to_plot, levels=levels, colors=colors)
            axs[i].set_xlabel(r'$\gamma_0$')
            axs[i].set_xticks([0.01, 0.5, 1])
            axs[i].set_xticklabels([0.01, 0.5, 1])
            axs[i].set_title(labels_[i])
            axs[i].set_xlim(0.01, 1)
            axs[i].set_ylim(0.01, 1)
        else:
            axs.set_aspect('equal')
            im = axs.tricontourf(triang, to_plot, levels=levels, colors=colors)
            axs.set_xlabel(r'$\gamma_0$')
            axs.set_xticks([0.01, 0.5, 1])
            axs.set_xticklabels([0.01, 0.5, 1])
            axs.set_title(labels_[i])
            axs.set_xlim(0.01, 1)
            axs.set_ylim(0.01, 1)

    if N_alphamodes>1:
        axs[0].set_ylabel(r'$S_0$')
        axs[0].set_yticks([0.01, 0.5, 1])
        axs[0].set_yticklabels([0.01, 0.5, 1])

        fig.subplots_adjust(bottom=0.2, top=0.95)
    else:
        axs.set_ylabel(r'$S_0$')
        axs.set_yticks([0.01, 0.5, 1])
        axs.set_yticklabels([0.01, 0.5, 1])


    return fig, axs, im, levels

def add_colorbar_to_plot_levels(fig, im, levels, ticks):
    if N_alphamodes > 1:
        cbar_ax = fig.add_axes([0.125, 0.15, 0.75, 0.02])
        cbar=fig.colorbar(im, cax=cbar_ax, orientation='horizontal', format='%.2e')
    else:
        cbar=fig.colorbar(im, orientation='horizontal', format='%.2e', aspect=50)
    cbar.set_ticks([a +0.5 for a in levels])
    cbar.set_ticklabels(ticks)
    return cbar


def closest_element_in_list(el, liste):
    closest_index=0
    distance=abs(el-liste[0])
    for i in range(1, len(liste)):
        dist = abs(el-liste[i])
        if(dist<=distance):
            distance=dist
            closest_index=i
    return liste[closest_index]

def index_closest_element_in_list(el, liste):
    closest_index=0
    distance=abs(el-liste[0])
    for i in range(1, len(liste)):
        dist = abs(el-liste[i])
        if(dist<=distance):
            distance=dist
            closest_index=i
    return closest_index

def remove_strings_from_file(matrices_folder,filename):
    file = open(filename + '.out', "r")
    metadata = []
    for x in file:
        name = x.replace(matrices_folder + '/RandTrix_Nr', '')
        name = name.replace('_Nc', ' ')
        name = name.replace('_Nest', ' ')
        name = name.replace('_Conn', ' ')
        name = name.replace('.txt', '')
        metadata.append(name)
    file.close()
    f = open(filename + '_filtered.out', 'w')
    for a in metadata:  # python will convert \n to os.linesep
        f.write(a + '\n')
    f.close()

# the idea is to plot points as full surface of a given colour on an axis
def plot_points_as_area(ax, points, col):
    points=np.array(points)
    popt, pcov, g0fit, S0fit = get_fit_cfr(points)
    ax.fill_between(g0fit,0, np.multiply(S0fit, g0fit), facecolor=col)
    return

def linear_function(x, a, b):
    return a*x+b

def exponential_function(x, a, b, c):
    return a*np.exp(-b*x)-c

def power_function(x, a, b, c, d):
    return a*np.power(x-d,b)+c

def fit_data(function_to_use, xdata, ydata):
    Npoints = len(xdata)
    if function_to_use==exponential_function:
        bounds=((0, 0, 0), (np.inf, np.inf, np.inf))
        popt, pcov = curve_fit(function_to_use, xdata, ydata, bounds=bounds, maxfev=9000000)
    elif function_to_use==power_function:
        bounds=((-np.inf, -np.inf, -np.inf, -np.inf), (np.inf, np.inf, np.inf, np.inf))
        popt, pcov = curve_fit(function_to_use, xdata, ydata, bounds=bounds, maxfev=9000000)
    elif function_to_use==linear_function:
        if Npoints>2:
            popt, pcov = curve_fit(function_to_use, xdata, ydata)
        else:
            a_ = (ydata[0]-ydata[1])/(xdata[0]-xdata[1])
            b_ = ydata[0]-a_*xdata[0]
            popt=[a_,b_]
            pcov=[[0,0],[0,0]]
    else:
        popt, pcov = curve_fit(function_to_use, xdata, ydata)
    fitted_y = [function_to_use(x, *popt) for x in xdata]
    return fitted_y, popt, np.sqrt(np.diagonal(pcov))


def zero_from_fit(function_to_use, popt, pcov):
    result = 0
    rel_error = 0
    error = 0

    if function_to_use==linear_function:
        rel_error = (pcov[0]/popt[0]+pcov[1]/popt[1])
        result = -popt[1]/popt[0]
    elif function_to_use==exponential_function:
        rel_error = pcov[0]/popt[0]+pcov[1]/popt[1]+pcov[2]/popt[2]
        result = np.log(abs(popt[0])/abs(popt[2]))/abs(popt[1])
    elif function_to_use==power_function:
        a,b,c,d = popt
        result = np.power(-c/a, 1./b)+d
    error = abs(rel_error*result)
    return result, error


def data_levels(data):
    NR=data[:,0]
    NS=data[:,1]
    nestedness=data[:,2]
    connectance=data[:,3]
    gamma0=data[:, 4::3]
    S0=data[:, 5::3]
    feasability=data[:, 6::3]
    ff_indices=[]
    for i in range(len(feasability)):
        # find the fully feasible_indices for this matrix at this alpha0
        full_feas_indices=[j for j in range(len(feasability[i])) if feasability[i,j]==1.]
        ff_indices.append(full_feas_indices)
    feasibility_level=[]
    # a level of -1 means you are not feasible, level 0 means you are feasible
    # for alpha0[0] but not after, level[1] means you are feasible for alpha0[1] but
    # not after and so on
    for i in range(len(feasability[0])):
        level=-1
        j=0
        exit=False
        while(not(exit) and j < len(feasability)):
            if i in ff_indices[j]:
                level+=1
            else:
                exit=True
            j+=1
        feasibility_level.append(level)
    return feasibility_level


# the first dimension are the different alpha0's
def plot_region_for_one_matrix(fig, ax, data):
    feasibility_level=feasibility_levels(data)
    max_level=max(feasibility_level)
    levels=[i for i in range(1,max_level+2)]

    # isbad=np.less(feasibility_level, 0.5)
    triang = tr.Triangulation(gamma0[0], S0[0])
    # mask = np.all(np.where(isbad[triang.triangles], True, False), axis=1)
    # triang.set_mask(mask)
    im=ax.tricontourf(triang, feasibility_level,levels=levels, cmap='jet_r')
    #im=ax.scatter(gamma0[0], S0[0], c=feasibility_level)
    #ax.set_yticklabels([])

    ax.set_xlim(0., 1.)
    ax.set_ylim(0., 1.)
    ax.set_xlabel(r'$\gamma_0$')
    ax.set_xticks([0, 0.5, 1])
    ax.set_xticklabels([0, 0.5, 1])
    ax.set_yticks([])
    #cbar.set_ticklabels(alpha0[0:max_level])
    return im, max_level

# data : get
def shrink_volume_for_one_matrix(data):
    NR=data[:,0]
    NS=data[:,1]
    nestedness=data[:,2]
    connectance=data[:,3]
    gamma0=data[:, 4::3]
    S0=data[:, 5::3]
    feasability=data[:, 6::3]
    ff_indices=[]
    volume=[]
    Nalpha0=len(feasability)
    Npoints=len(feasability[0])

    for i in range(Nalpha0):
        # find the fully feasible_indices for this matrix at this alpha0
        local_volume=0.
        for el in feasability[i]:
            local_volume+=el
        #full_feas_indices=[j for j in range(len(feasability[i])) if feasability[i,j]==1.]
        volume.append(local_volume/Npoints)

    return volume

# scurve contains the data of the shrinkage curve
def fit_shrinkage_curve(alpha0_,scurve_, fit_function_, take_zero_points_):
    # we do a linear fit removing the first point up ontil the first zero point
    end_index=len(scurve_)
    zero_indices = [j for j in range(len(scurve_)) if scurve_[j]==0.]
    start_index=0
    if(len(zero_indices)>0):
        end_index=np.min(zero_indices)
    if take_zero_points_:
        end_index=len(scurve_)
    if end_index>2:
        fitted_curve, popt, perr = fit_data(fit_function_, alpha0_[start_index:end_index], scurve_[start_index:end_index])
        if fit_function_==linear_function:
            from scipy.stats import linregress
            slope, intercept, r_value, p_value, stderr = linregress(alpha0_[start_index:end_index], scurve_[start_index:end_index])
            print("p-value: ", p_value)
            estimated_decay_rate=-slope
        if fit_function_==exponential_function:
            estimated_decay_rate=(popt[1], perr[1])

        if fit_function_==power_function:
            estimated_decay_rate=popt[1]
        # we now get the estimated critical alpha0 and the estimated vol at alpha0=0
        estimated_alpha_crit, err_alpha_crit=zero_from_fit(fit_function_, popt, perr)
        estimated_vol_zero_syntrophy = fit_function_(0, *popt)
    else:
        estimated_vol_zero_syntrophy = scurve_[0]
        fitted_curve=scurve_[1:end_index]
        if end_index==2:
            a = (scurve_[0]-scurve_[1])/(alpha0_[0]-alpha0_[1])
            b = scurve_[0]-a*alpha0_[0]
            estimated_alpha_crit= -b/a
            estimated_decay_rate=-a
        else:
            estimated_alpha_crit=alpha0_[1]
            estimated_decay_rate=-(scurve_[1]-scurve_[0])/(alpha0_[1]-alpha0_[0])
    remaining_points=[0 for j in range(end_index, len(scurve_))]
    if start_index==1:
        fitted_curve = [estimated_vol_zero_syntrophy]+fitted_curve
    fitted_alpha0 = alpha0_[start_index:end_index]
    print(fitted_alpha0, fitted_curve)
    return fitted_alpha0, fitted_curve, estimated_alpha_crit, estimated_vol_zero_syntrophy, estimated_decay_rate



def func_to_fit(x, k2, k3):
    return k2*x+k3

def get_fit_cfr(points):
    gamma0=points[:,0]
    S0=points[:,1]

    # for a given gamma0 we find the largest S0 possible -> part of the border
    border_points=[]
    g0=sorted(list(set(gamma0)))
    for g in g0:
        S = np.array([S0[i] for i in range(len(points)) if gamma0[i]==g])
        border_points.append([g, np.max(S)])
    S = sorted(list(set(S0)))
    border_points=np.array(border_points)
    gborder=border_points[:,0]
    Sborder=border_points[:,1]

    popt, pcov = curve_fit(func_to_fit, Sborder, np.multiply(gborder, Sborder))
    fitted_g0S0 = [func_to_fit(a, *popt) for a in np.multiply(gborder, Sborder)]

    fitted_S0 = np.array([fitted_g0S0[i]/gborder[i] for i in range(len(gborder))])

    return popt, pcov, gborder, fitted_S0

# finds the asymptote value of an asymptotic function y, which is assumed to be constant
# or constantly oscillating for its last Npoints
def asymptote(y, Npoints):
    return np.mean(y[-Npoints:])
def remove_strings_from_file(matrices_folder,filename):
    file = open(filename + '.out', "r")
    metadata = []
    for x in file:
        name = x.replace(matrices_folder + '/RandTrix_Nr', '')
        name = name.replace('_Nc', ' ')
        name = name.replace('_Nest', ' ')
        name = name.replace('_Conn', ' ')
        name = name.replace('.txt', '')
        metadata.append(name)
    file.close()
    f = open(filename + '_filtered.out', 'w')
    for a in metadata:  # python will convert \n to os.linesep
        f.write(a + '\n')
    f.close()
    return

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

def extract_numbers_from_string(filename):
    # printing original string
    # using re.findall()
    # getting numbers from string
    temp = re.findall(r"[-+]?\d*\.\d+|\d+", filename)
    res = list(map(float, temp))
    return res

def mat_connectance(mat):
    N = len(mat)
    M = len(mat[0])
    links=0

    for i in range(N):
        for j in range(M):
            if mat[i][j]!=0:
                links+=1.
    return links/(N*M)

def mat_nestedness(mat):
    M = []

    n_rows = len(mat)
    n_cols = len(mat[0])

    for i in range(n_rows):
        M.append([])
        for j in range(n_cols):
            if mat[i][j]!=0:
                M[i].append(1)
            else:
                M[i].append(0)

    n = []
    for i in range(n_rows):
        n.append(0.)
        for k in range(n_cols):
            n[i] += M[i][k]

    n_tilde = []
    for i in range(n_rows):
        n_tilde.append([])
        for j in range(n_cols):
            n_tilde[i].append(0.)
            for k in range(n_cols):
                n_tilde[i][j]+=M[i][k]*M[j][k]

    nest = 0.
    num = sum([n_tilde[i][j] for j in range(n_cols) for i in range(n_rows) if i < j])
    denom = sum([min([n[i], n[j]]) for j in range(n_cols) for i in range(n_rows) if i < j])

    if denom==0:
        nest=0
    else:
        nest = num/denom


    return nest

def eco_energy(A, G, alpha0=1, gamma0=1, R0=1):
    NR=len(A)
    AG=A@G
    GG=np.transpose(G)@G
    Z =np.zeros(NR)

    energy=0.


    for mu in range(NR):
        Z[mu]+=(alpha0*AG[mu][mu]-gamma0*R0*GG[mu][mu]);
        for nu in range(NR):
            if nu!=mu:
                Z[mu]+=abs(alpha0*AG[mu][nu]-gamma0*R0*GG[mu][nu])



    coeff = 1;
    for mu in range(NR):
        energy+=coeff*Z[mu]
    return energy;

def compute_feasibility_data(to_compute, volume_data, decay_rate_data):
    return compute_data(to_compute, volume_data, decay_rate_data, data_type="feasible")

def compute_lds_data(to_compute, volume_data, decay_rate_data):
    return compute_data(to_compute, volume_data, decay_rate_data, data_type="ld stable")

def compute_largest_eigenvalue_data(alpha_mode_, alpha0_, filename, optimal_LRI_folder, consumption_matrix_folder, to_compute, save_file):
    columns = ["NR", "NS", "connG", "nestG", "alpha_mode"]
    filter_data(alpha_mode_,alpha0_,filename,optimal_LRI_folder, consumption_matrix_folder)
    data_region = load_data_region(alpha_mode_, alpha0_, filename, optimal_LRI_folder, type_=np.complex)
    if 'av. dominant eigenvalue' in to_compute:
        df = pd.DataFrame(columns=columns+["alpha0", "av. dominant eigenvalue"])
        for amode in range(len(alpha_mode_)):
            for a0 in range(len(alpha0_)):
                for mat in range(len(data_region[amode, a0])):
                    data = np.real(np.ma.masked_invalid(data_region[amode, a0, mat, 6::3]))
                    av_dom_eig = np.NaN
                    if len(data[~data.mask]) > 0:
                        av_dom_eig = np.sum(data[~data.mask])/len(data[~data.mask])
                    dict={
                        "NR": int(data_region[amode, a0, mat,0]),
                        "NS": int(data_region[amode, a0, mat,1]),
                        "connG": closest_element_in_list(data_region[amode, a0, mat,3], all_connectance),
                        "nestG": closest_element_in_list(data_region[amode, a0, mat,2], all_nestedness),
                        "alpha_mode":alpha_mode_[amode],
                        "alpha0": alpha0_[a0],
                        "av. dominant eigenvalue": np.abs(av_dom_eig)
                    }
                    df=df.append(pd.DataFrame([dict]), sort=False)
        df.to_csv(save_file, index=False)
        print("Saved av. dominant eigenvalues in file ", save_file)
    return

def compute_data(to_compute, volume_data, decay_rate_data, data_type):
    columns = ["NR", "NS", "connG", "nestG", "alpha_mode"]

    if data_type+' volume' in to_compute:

        alpha_mode = volume_data['alpha_mode']
        alpha0 = volume_data['alpha0']
        filename = volume_data['data file']
        optimal_LRI_folder = volume_data['OM folder']
        consumption_matrix_folder = volume_data['G matrices folder']
        save_file = volume_data['save file']

        filter_data(alpha_mode, alpha0, filename, optimal_LRI_folder, consumption_matrix_folder)
        data_region = load_data_region(alpha_mode, alpha0, filename, optimal_LRI_folder)
        alpha0=np.array(alpha0)

        df = pd.DataFrame(columns=columns+["alpha0", data_type+" volume"])
        for al_mode in range(len(alpha_mode)):
            for a0 in range(len(alpha0)):
                for mat in range(len(data_region[al_mode, a0])):
                    data=data_region[al_mode, a0, mat, 6::3]
                    volume=np.sum(data)/len(data)
                    dict={
                        "NR": int(data_region[al_mode, a0, mat,0]),
                        "NS": int(data_region[al_mode, a0, mat,1]),
                        "connG": closest_element_in_list(data_region[al_mode, a0, mat,3], all_connectance),
                        "nestG": closest_element_in_list(data_region[al_mode, a0, mat,2], all_nestedness),
                        "alpha_mode":alpha_mode[al_mode],
                        "alpha0": alpha0[a0],
                        data_type+" volume": volume,
                    }

                    df=df.append(pd.DataFrame([dict]), sort=False)
        df.to_csv(save_file, index=False)
        print("Saved "+data_type+" volumes in file ", save_file)

    if data_type+' decay rate' in to_compute:

        alpha_mode = decay_rate_data['alpha_mode']
        alpha0 = decay_rate_data['alpha0_range']
        filename = decay_rate_data['data file']
        optimal_LRI_folder = decay_rate_data['OM folder']
        consumption_matrix_folder = decay_rate_data['G matrices folder']
        save_file = decay_rate_data['save file']

        filter_data(alpha_mode, alpha0, filename, optimal_LRI_folder, consumption_matrix_folder)
        data_region = load_data_region(alpha_mode, alpha0, filename, optimal_LRI_folder)

        df = pd.DataFrame(columns=columns+[data_type+" decay rate", data_type+" decay rate error"])
        for al_mode in range(len(alpha_mode)):
            data_vols=[]
            for a0 in range(len(alpha0)):
                local_vols = []
                for mat in range(len(data_region[al_mode, a0])):
                    data=data_region[al_mode, a0, mat, 6::3]
                    data_volume=np.sum(data)/len(data)
                    local_vols.append(data_volume)
                data_vols.append(local_vols)
            data_vols=np.transpose(data_vols)
            for mat in range(len(data_vols)):
                dict={
                    "NR": int(data_region[al_mode, 0, mat,0]),
                    "NS": int(data_region[al_mode, 0, mat,1]),
                    "connG": closest_element_in_list(data_region[al_mode, 0, mat,3], all_connectance),
                    "nestG": closest_element_in_list(data_region[al_mode, 0, mat,2], all_nestedness),
                    "alpha_mode": alpha_mode[al_mode]
                    }
                fitted_y, popt, perr = fit_data(exponential_function, alpha0, data_vols[mat])
                dict[data_type+' decay rate']=popt[1]
                dict[data_type+' decay rate error']=perr[1]
                df=df.append(pd.DataFrame([dict]), sort=False)
        df.to_csv(save_file, index=False)
        print("Saved "+data_type+" decay rates in file ", save_file)
    return

def plot_feasible_volume(ax, data_file, width, shift, alpha0_, alpha_mode_):
    ax = plot_volumes(ax, data_file, width, shift, alpha0_, alpha_mode_, data_type='feasible')
    ax.set_ylabel('Feasible volume')
    return ax
def plot_lds_volume(ax, data_file, width, shift, alpha0_, alpha_mode_):
    ax = plot_volumes(ax, data_file, width, shift, alpha0_, alpha_mode_, data_type='ld stable')
    ax.set_ylabel('Dynamically stable volume')
    return ax
def plot_largest_eigenvalue(ax, data_file, width, shift, alpha0_, alpha_mode_):
    ax = plot_volumes(ax, data_file, width, shift, alpha0_, alpha_mode_, data_type='av. dominant eigenvalue')
    ax.set_ylabel('-Dominant eigenvalue')
    return ax
def plot_feasible_decay_rates(ax, decay_rate_data):
    return plot_decay_rates(ax, decay_rate_data, type='feasible')

def plot_p_values_vs_alpha0(axes, volume, type):
    data = pd.read_csv(volume['data file'])
    # take all possibles alpha_modes pairs
    a0modes = volume['alpha mode']
    pairs = [(a0modes[i], a0modes[j]) for i in range(len(a0modes)) for j in range(i+1, len(a0modes))]
    for pair in pairs:
        a0mode1 = pair[0]
        a0mode2 = pair[1]
        p_vals = []
        a_vals = []
        for a0 in volume['alpha0']:
            indices1 = [x for x in range(data.shape[0]) if (data.loc[x]['alpha_mode']==a0mode1 and closest_element_in_list(data.loc[x]['alpha0'], volume['alpha0'])==a0)]
            indices2 = [x for x in range(data.shape[0]) if (data.loc[x]['alpha_mode']==a0mode2 and closest_element_in_list(data.loc[x]['alpha0'], volume['alpha0'])==a0)]

            distrib1 = data.loc[indices1][type+' volume'].to_numpy()
            distrib2 = data.loc[indices2][type+' volume'].to_numpy()

            statistics, pval = wilcoxon(x=distrib1, y=distrib2, alternative='less', zero_method='zsplit')
            p_vals.append(pval)
            a_vals.append(a0)
        axes.plot(a_vals, p_vals, label=alpha_mode_label[a0mode1]+'-'+alpha_mode_label[a0mode2])
        axes.legend()
        axes.set_xlabel(r'$\alpha_0$')
        axes.set_ylabel(r'$p$-value ('+type+')')
        axes.set_title(r'Null hypothesis for pair P1-P2 : median(P1) $>$ median(P2). Small $p$ : null hyp. should be rejected', fontsize=10)
    return axes

def plot_data(figures_to_plot, volume, decay_rates, type):
    if type+' volume' in figures_to_plot:
        fig = plt.figure(type+' volume')
        ax = fig.add_subplot(111)
        intrashift = volume['intrashift']
        intershift = volume['intershift']
        shift = [intershift, intrashift]
        width=volume['width']
        if type=='feasible':
            ax = plot_feasible_volume(ax, volume['data file'], width, shift,
                        volume['alpha0'], volume['alpha mode'])
        elif type=='ld stable':
            ax = plot_lds_volume(ax, volume['data file'], width, shift,
                        volume['alpha0'], volume['alpha mode'])
        fig.tight_layout()
        fig.savefig(volume['save name'])

    if type+' decay rate' in figures_to_plot:
        Namodes = len(decay_rates['alpha mode'])
        axes = []
        figs = []
        for i in range(Namodes):
            figs.append(plt.figure())
            figs.append(plt.figure())
            axes.append([figs[2*i].add_subplot(111), figs[2*i+1].add_subplot(111)])
        axes = plot_decay_rates(axes, decay_rates, type)
        for i in range(Namodes):
            figs[2*i].tight_layout()
            figs[2*i+1].tight_layout()

            figs[2*i].savefig(decay_rates['save name']+'_'+decay_rates['alpha mode'][i]+'_fixed_conn.pdf')
            figs[2*i+1].savefig(decay_rates['save name']+'_'+decay_rates['alpha mode'][i]+'_fixed_nest.pdf')

    if type+' p-value' in figures_to_plot:
        Namodes = len(decay_rates['alpha mode'])
        fig = plt.figure()
        axes = fig.add_subplot(111)
        axes = plot_p_values_vs_alpha0(axes, volume, type)
        fig.tight_layout()
        fig.savefig(decay_rates['save name']+'_p_value.pdf')
    return


def plot_decay_rates(axs, decay_rate_data, type):
    data = pd.read_csv(decay_rate_data['data file'])
    for j in range(len(decay_rate_data['alpha mode'])):
        amode = decay_rate_data['alpha mode'][j]
        for i in range(len(all_connectance)):
            connG = all_connectance[i]
            indices = [x for x in range(data.shape[0]) if closest_element_in_list(data.loc[x]['connG'], all_connectance)==connG and data.loc[x]['alpha_mode']==amode]
            to_plot = data.loc[indices]
            to_plot = to_plot.sort_values(by='nestG')
            to_plot.plot(x = 'nestG', y=type+' decay rate', yerr=type+' decay rate error', ax=axs[j][1], label=r'$\kappa_G \approx'+str(connG)+'$',
                    linestyle='solid', color=conn_colours[i], marker=alpha_mode_sym[amode])
        min_nest = np.min(data['nestG'])
        max_nest = np.max(data['nestG'])
        axs[j][1].set_xlim(min_nest-0.1*(max_nest-min_nest), max_nest+0.1*(max_nest-min_nest))
        axs[j][1].set_xlabel(labels['nestG'])
        axs[j][1].set_ylabel(labels[type+' decay rate'])
        axs[j][1].legend(bbox_to_anchor=(1.05,1.), loc='upper left')
        axs[j][1].set_title(alpha_mode_label[amode])

    for j in range(len(decay_rate_data['alpha mode'])):
        amode = decay_rate_data['alpha mode'][j]
        for i in range(len(all_nestedness)):
            nestG = all_nestedness[i]
            indices = [x for x in range(data.shape[0]) if closest_element_in_list(data.loc[x]['nestG'], all_nestedness)==nestG and data.loc[x]['alpha_mode']==amode]
            to_plot = data.loc[indices]
            to_plot = to_plot.sort_values(by='connG')
            to_plot.plot(x = 'connG', y=type+' decay rate', yerr=type+' decay rate error', ax=axs[j][0], label=r'$\eta_G \approx'+str(nestG)+'$',
                    linestyle='solid', color=nest_colours[i], marker=alpha_mode_sym[amode])
        min_conn = np.min(data['connG'])
        max_conn = np.max(data['connG'])
        axs[j][0].set_xlim(min_conn-0.1*(max_conn-min_conn), max_conn+0.1*(max_conn-min_conn))
        axs[j][0].set_xlabel(labels['connG'])
        axs[j][0].set_ylabel(labels[type+' decay rate'])
        axs[j][0].legend(bbox_to_anchor=(1.05,1.), loc='upper left')
        axs[j][0].set_title(alpha_mode_label[amode])

    return axs

def plot_volumes(ax, data_file, width, shift, alpha0_, alpha_mode_, data_type):
    df = pd.read_csv(data_file)

    intershift = shift[0]
    intrashift = shift[1]

    N_alphamodes = len(alpha_mode_)
    N_alpha0 = len(alpha0_)
    L = N_alphamodes*(width+intrashift)-intrashift+intershift
    xs = 0.5*(width+intershift)

    for a in df['alpha0']:
        df['alpha0'] = df['alpha0'].replace([a], closest_element_in_list(a, alpha0_))

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
        data=df[df['alpha_mode']==amode]
        to_plot=pd.DataFrame(columns=['matrix'])
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
        for a0 in alpha0_:
            string = data_type+' volume'
            if data_type == "av. dominant eigenvalue":
                string = 'av. dominant eigenvalue'
            volumes = data[data['alpha0']==a0][string].to_numpy()
            volumes = volumes[~np.isnan(volumes)]
            means.append(np.mean(volumes))
            to_plot=to_plot.append(pd.DataFrame([volumes]), sort=False)
        all_means.append(means)
        positions = np.linspace(start=xs+i*(width+intrashift), stop=xs+(N_alpha0-1)*L+i*(width+intrashift), num=N_alpha0)
        all_positions.append(positions)
        ax.boxplot(to_plot, positions=positions, widths=width,
                patch_artist=True , boxprops=boxprops,
                meanprops=meanprops, medianprops=medianprops,
                flierprops=flierprops, whiskerprops=whiskerprops,
                capprops=capprops)
        legend_els.append(Patch(facecolor=facecol, edgecolor=facecol, linewidth=0, label=alpha_mode_label[amode]))

    for i in range(N_alphamodes):
        amode = alpha_mode_[i]
        col = alpha_mode_colours[amode]
        marker = alpha_mode_sym[amode]
        ax.plot(all_positions[i], all_means[i], marker=marker, markerfacecolor='black', markersize=5,
                linestyle='solid', color=col, markeredgecolor='black', markeredgewidth=0.5)



    ax.set_xticks(ticks)
    ax.set_xticklabels([a*1e3 for a in alpha0_], rotation=45)
    ax.set_xlim(0, xs+(N_alpha0-1)*L+(N_alphamodes-1)*(width+intrashift)+(width+intershift)*0.5)
    ax.legend(handles=legend_els, bbox_to_anchor=(0., 1.02, 1., .102), loc='center',
           ncol=len(legend_els), borderaxespad=0., fontsize=12)
    ax.set_title('')
    ax.set_yscale('linear')
    ax.set_xlabel(r'$\alpha_0 \times 10^{3}$')
    return ax
