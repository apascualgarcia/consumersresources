import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits import mplot3d
import matplotlib
# rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
# # for Palatino and other serif fonts use:
# # rc('font',**{'family':'serif','serif':['Palatino']})
# rc('text', usetex=False)

folder = './data_output'
filename = 'cd_alpha_max_configuration_comparison_NR25_NS25_s05_a1_critical_matrix_list_NR25_NS25_3'
save_folder = './plots'
save_name = filename
title = r''
errorbar = 'no_errorbar'
matrices_folder = './matrices/Nr25_Nc25'
actual_connectance = [0.08, 0.13, 0.18, 0.23, 0.28, 0.33, 0.38, 0.43]

def remove_strings_from_file(filename):
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

def nearest_el_in_list(x, liste):
    result = liste[0]
    for i in range(len(liste)):
        if((x-liste[i])**2 < (x-result)**2):
            result = liste[i]
    return result


remove_strings_from_file(folder + '/' + filename)
# each line of the data corresponds to the points related to a specific matrix
data = np.loadtxt(folder + '/' + filename + '_filtered.out')
# metadata = data[:,:4]
# for j in range(len(data)):
#     # we round the nestedness to two decimal values
#     data[j,2] = round(data[j,2], 2)
#     # we filter the connectance such that it's rounded to the actual target value
#     data[j,3] = nearest_el_in_list(data[j,3], actual_connectance)

metadata = data[:4]
data[2] = round(data[2],2)
data[3] = nearest_el_in_list(data[3], actual_connectance)

alpha = np.array([data[i] for i in range(4,len(data),3)])
delta_crit=np.array([data[i] for i in range(5, len(data),3)])
delta_error=np.array([data[i] for i in range(6, len(data),3)])

fig1 = plt.figure('Critical delta varying syntrophy')
ax1 = fig1.add_subplot(111)
ax1.errorbar(x=alpha, y=delta_crit, yerr=delta_error,fmt='+-', linestyle='dotted', elinewidth=1)
ax1.set_title(r'$\sigma_0=0.5$, $N_R=N_S=25$, $\kappa=$'+str(data[3])+', $\eta=$'+str(data[2]))
ax1.set_xlabel(r'$\alpha_0$')
ax1.set_ylabel(r'$\Delta^*$')
fig1.tight_layout()
fig1.savefig(save_folder+'/'+filename)
