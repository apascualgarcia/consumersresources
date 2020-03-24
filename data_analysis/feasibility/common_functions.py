import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.tri as tr

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

# the first dimension are the different alpha0's
def plot_feasibility_region_for_one_matrix(data):
    NR=data[:,0]
    NS=data[:,1]
    nestedness=data[:,2]
    connectance=data[:,3]
    gamma0=data[:, 4::3]
    S0=data[:, 5::3]
    feasability=data[:, 6::3]
    ff_indices=[]
    fig=plt.figure()
    ax=fig.add_subplot(111)
    for i in range(len(feasability)):
        # find the fully feasible_indices for this matrix at this alpha0
        full_feas_indices=[j for j in range(len(feasability[i])) if feasability[i,j]==1.]
        ff_indices.append(full_feas_indices)
    feasibility_level=[]
    for i in range(len(feasability[0])):
        level=0
        j=0
        exit=False
        while(not(exit) and j < len(feasability)):
            if i in ff_indices[j]:
                level+=1
            else:
                exit=True
            j+=1
        feasibility_level.append(level)
    max_level=max(feasibility_level)
    levels=[i for i in range(max_level+1)]

    isbad=np.less(feasibility_level, 0.5)
    triang = tr.Triangulation(gamma0[0], S0[0])
    mask = np.all(np.where(isbad[triang.triangles], True, False), axis=1)
    triang.set_mask(mask)
    im=ax.tricontourf(triang, feasibility_level, level=levels, cmap='Greys')
    #im=ax.scatter(gamma0[0], S0[0], c=feasibility_level)

    ax.set_xlim(0., 1.)
    ax.set_ylim(0., 1.)
    ax.set_xlabel(r'$\gamma_0$')
    ax.set_ylabel(r'$S_0$')
    cbar=fig.colorbar(im)
    cbar.set_ticks(levels[0:max_level])
    #cbar.set_ticklabels(alpha0[0:max_level])
    ax.set_title(r'$N_R=25, N_S=25$ ')
    fig.tight_layout()
    #plt.show()
    return



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
