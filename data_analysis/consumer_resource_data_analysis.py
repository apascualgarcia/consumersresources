import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.tri as tr
import matplotlib as mpl
import copy
import sys
np.set_printoptions(threshold=sys.maxsize)
mpl.rcParams['lines.linewidth']=2.5
mpl.rcParams['lines.markersize']=10
mpl.rcParams['lines.markeredgewidth']=3

alpha_mode=['fully_connected', 'no_release_when_eat', 'optimal_matrix', 'random_structure']
alpha_mode_colours=['blue', 'green', 'red', 'black']
label=['fully connected', 'no intraspecific syntrophy', 'LRI regime', 'random structure']
alpha0=[0, 1.3e-3, 2.6e-3, 3.9e-3, 5.2e-3, 6.5e-3, 7.8e-3, 9.1e-3, 1.04e-2, 1.4e-2]
all_nestedness=[0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]
all_connectance=[0.08, 0.13, 0.18, 0.23, 0.28, 0.33, 0.38, 0.43]
min_gamma0, max_gamma0=0.01, 1
min_S0, max_S0=0.01, 1
nestedness_label=r'Ecological overlap $\eta_G$'
connectance_label=r'Connectance $\kappa_G$'

N_alphamodes=len(alpha_mode)

# region contains all data for each alphamode, alpha0 and network
def plot_common_region(region, alpha_mode_, colors_, labels_):
    fig, axs = plt.subplots(1, N_alphamodes, sharey=True, sharex=True, figsize=(3.5*N_alphamodes,4.5))
    quantity=region[:,:,:,6::3]
    gamma0=region[:,:,:,4::3]
    S0=region[:,:,:,5::3]
    Nmatrices=len(region[0,0])

    # do that computation for each alpha_mode
    for k in range(len(region)):
        ax=axs[k]

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
    axs[0].set_ylabel(r'$S_0$')
    axs[0].set_yticks([0, 0.5, 1])
    axs[0].set_yticklabels([0, 0.5, 1])

    fig.subplots_adjust(bottom=0.2, top=0.95)

    return fig, axs, im, levels

def filter_data(alpha_mode_, alpha0_, filename_, optimal_LRI_folder_, consumption_matrix_folder_):
    for al_mo in alpha_mode_:
        for a in alpha0_:
            file=filename_+'_'+al_mo+'_'+optimal_LRI_folder_+'_alpha0='+str(a)
            remove_strings_from_file(consumption_matrix_folder_, file)
    return
def load_data_region(alpha_mode_, alpha0_, filename_, optimal_LRI_folder_):
    region=[]
    for al_mo in alpha_mode_:
        local_vector=[]
        for a in alpha0_:
            file=filename_+'_'+al_mo+'_'+optimal_LRI_folder_+'_alpha0='+str(a)+'_filtered.out'
            local_data=np.loadtxt(file)
            print('Loading file', file, 'which contains the data of', len(local_data), 'matrices')
            local_vector.append(local_data)
        region.append(local_vector)
    region=np.array(region)
    print('Dimensions of tableau : ', len(region), 'x', len(region[0]),'x',len(region[0,0]))
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


    fig, axs = plt.subplots(1, N_alphamodes, sharey=True, sharex=True, figsize=(3.5*N_alphamodes,4.5))
    function_levels=levels_different_alpha0(data)
    max_level=np.amax(function_levels)
    levels=[i for i in range(max_level+2)]
    triang = tr.Triangulation(gamma0, S0)
    for i in range(len(data)):
        to_plot = function_levels[i]
        axs[i].set_aspect('equal')
        im = axs[i].tricontourf(triang, to_plot, levels=levels, colors=colors)
        axs[i].set_xlabel(r'$\gamma_0$')
        axs[i].set_xticks([0.01, 0.5, 1])
        axs[i].set_xticklabels([0.01, 0.5, 1])
        axs[i].set_title(labels_[i])
        axs[i].set_xlim(0.01, 1)
        axs[i].set_ylim(0.01, 1)
    axs[0].set_ylabel(r'$S_0$')
    axs[0].set_yticks([0.01, 0.5, 1])
    axs[0].set_yticklabels([0.01, 0.5, 1])

    fig.subplots_adjust(bottom=0.2, top=0.95)

    return fig, axs, im, levels

def add_colorbar_to_plot_levels(fig, im, levels, ticks):
    cbar_ax = fig.add_axes([0.125, 0.15, 0.75, 0.02])
    cbar=fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
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

def power_function(x, a, b, c):
    return a*np.power(x,b)+c

def fit_data(function_to_use, xdata, ydata):
    Npoints = len(xdata)
    if function_to_use==exponential_function:
        bounds=((0, 0, 0), (np.inf, np.inf, np.inf))
        popt, pcov = curve_fit(function_to_use, xdata, ydata, bounds=bounds, maxfev=9000000)
    elif function_to_use==power_function:
        bounds=((-np.inf, -np.inf, -np.inf), (np.inf, np.inf, np.inf))
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
        a,b,c = popt
        result = np.power(-c/a, 1./b)
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
        full_feas_indices=[j for j in range(len(feasability[i])) if feasability[i,j]==1.]
        volume.append(len(full_feas_indices)/Npoints)

    return volume

# scurve contains the data of the shrinkage curve
def fit_shrinkage_curve(alpha0_,scurve_):
    # we do a linear fit removing the first point up ontil the first zero point
    end_index=len(scurve_)
    zero_indices = [j for j in range(len(scurve_)) if scurve_[j]==0.]
    start_index=1
    if(len(zero_indices)>0):
        end_index=np.min(zero_indices)
    if end_index>2:
        fit_function=linear_function
        fitted_curve, popt, perr = fit_data(fit_function, alpha0_[start_index:end_index], scurve_[start_index:end_index])

        from scipy.stats import linregress
        slope, intercept, r_value, p_value, stderr = linregress(alpha0_[start_index:end_index], scurve_[start_index:end_index])
        print("p-value: ", p_value)
        # we now get the estimated critical alpha0 and the estimated vol at alpha0=0
        estimated_alpha_crit, err_alpha_crit=zero_from_fit(fit_function, popt, perr)
        estimated_vol_zero_syntrophy = fit_function(0, *popt)
    else:
        estimated_vol_zero_syntrophy = scurve_[0]
        fitted_curve=scurve_[1:end_index]
        if end_index==2:
            a = (scurve_[0]-scurve_[1])/(alpha0_[0]-alpha0_[1])
            b = scurve_[0]-a*alpha0_[0]
            estimated_alpha_crit= -b/a
        else:
            estimated_alpha_crit=alpha0
    remaining_points=[0 for j in range(end_index, len(scurve_))]
    fitted_curve = [estimated_vol_zero_syntrophy]+fitted_curve
    fitted_alpha0 = alpha0_[0:end_index]

    return fitted_alpha0, fitted_curve, estimated_alpha_crit, estimated_vol_zero_syntrophy



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
