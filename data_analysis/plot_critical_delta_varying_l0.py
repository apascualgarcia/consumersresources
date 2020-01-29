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
filename = 'critical_delta_varying_l0_alpha0=05'
save_folder = './plots'
save_name = filename
connectance=0.13
nestedness=0.3
alpha0=0.5


title = r'$\sigma_0=0.5$, $N_R=N_S=25$, $\alpha_0='+str(alpha0)+'$ $\kappa=$'+str(connectance)+', $\eta=$'+str(nestedness)


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


data = np.loadtxt(folder + '/' + filename + '.out')

l0 = np.array([data[i,0] for i in range(len(data))])
delta_crit=np.array([data[i,1] for i in range(len(data))])
delta_error=np.array([data[i,2] for i in range(len(data))])

fig1 = plt.figure('Critical delta varying external feeding')
ax1 = fig1.add_subplot(111)
ax1.errorbar(x=l0, y=delta_crit, yerr=delta_error,fmt='+-', linestyle='dotted', elinewidth=1)
ax1.set_title(title)
ax1.set_xlabel(r'$\l_0$')
ax1.set_ylabel(r'$\Delta^*$')
ax1.set_yscale('linear')
fig1.tight_layout()
fig1.savefig(save_folder+'/'+filename)
