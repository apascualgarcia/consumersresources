import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.tri as tr

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
    if function_to_use==exponential_function:
        bounds=((0, 0, 0), (np.inf, np.inf, np.inf))
        popt, pcov = curve_fit(function_to_use, xdata, ydata, bounds=bounds, maxfev=9000000)
    elif function_to_use==power_function:
        bounds=((-np.inf, -np.inf, -np.inf), (np.inf, np.inf, np.inf))
        popt, pcov = curve_fit(function_to_use, xdata, ydata, bounds=bounds, maxfev=9000000)

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
    for i in range(len(feasability)):
        # find the fully feasible_indices for this matrix at this alpha0
        full_feas_indices=[j for j in range(len(feasability[i])) if feasability[i,j]==1.]
        volume.append(len(full_feas_indices))
    #init_vol = volume[0]
    init_vol=900
    for i in range(len(feasability)):
        volume[i]=volume[i]/init_vol

    return volume



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
