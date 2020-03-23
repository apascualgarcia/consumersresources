import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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
