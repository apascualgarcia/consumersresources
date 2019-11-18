import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
# rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
# # for Palatino and other serif fonts use:
# # rc('font',**{'family':'serif','serif':['Palatino']})
# rc('text', usetex=False)

folder = './data_output'
filename = 'critical_delta_NR25_NS25'
save_folder = './plots'
save_name = filename
title = r'Nr25_Nc25_Nest0.3_Conn0.1296 (1000 runs per point)'
fs = 12
error_bar_width = 1.
cap_width = 1.5
markeredgewidth = 2
markersize = 6
errorbar = 'no_errorbar'
matrices_folder = './matrices/Nr25_Nc25'


def remove_strings_from_file(filename):
    if(not os.path.isfile(filename + '_filtered.out')):
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


filename = folder + '/' + filename
remove_strings_from_file(filename)
data = np.loadtxt(filename + '_filtered.out')
metadata = data[:, :4]


# sorted_size[i] contains the different system sizes
sorted_size = []
# sorted_mprop[i] contains all matrices of size sorted_size[i]
sorted_mprop = []

for i in range(0, len(metadata)):
    size_system = metadata[i, :2].tolist()
    if size_system in sorted_size:
        index = sorted_size.index(size_system)
        sorted_mprop[index].append(data[i, 2:].tolist())
    else:
        sorted_size.append(size_system)
        sorted_mprop.append([data[i, 2:].tolist()])
sorted_mprop = np.array(sorted_mprop)
sorted_size = np.array(sorted_size)
print("Data sorted.")

for i in range(len(sorted_size)):
    [NR, NS] = sorted_size[i]
    data = sorted_mprop[i]
    # we filter the data to round the nestedness up to 2 decimals
    data[:, 1] = np.array([round(data[i, 1], 2) for i in range(len(data))])
    # lists the different connectances available
    connectance = np.array(sorted(list(set(data[:, 0]))))
    # lists the different nestednesses available
    nestedness = np.array(sorted(list(set(data[:, 1]))))

    fig1 = plt.figure('Nestedness for fixed connectance')
    ax1 = fig1.add_subplot(111)
    for conn in connectance:
        nest = np.array([data[i, 1]
                         for i in range(len(data)) if data[i, 0] == conn])
        dcrit = np.array([data[i, 2]
                          for i in range(len(data)) if data[i, 0] == conn])
        indices = np.argsort(nest)
        sorted_dcrit = np.array([dcrit[i] for i in indices])
        sorted_nest = np.sort(nest)
        ax1.plot(sorted_nest, sorted_dcrit, '+-', linestyle='dotted',
                 label=r'Connectance $\approx$' + str(conn), markersize=markersize, markeredgewidth=markeredgewidth)
    ax1.set_title(r'$N_R=$' + str(NR) + r', $N_S=$' + str(NS))
    ax1.set_xlabel(r'Nestedness', fontsize=fs)
    ax1.set_ylabel(r'$\Delta^*$', fontsize=fs)
    ax1.legend()

    fig2 = plt.figure('Connectance for fixed nestedness')
    ax2 = fig2.add_subplot(111)
    for nest in nestedness:
        conn = np.array([data[i, 0]
                         for i in range(len(data)) if data[i, 1] == nest])
        dcrit = np.array([data[i, 2]
                          for i in range(len(data)) if data[i, 1] == nest])
        indices = np.argsort(conn)
        sorted_dcrit = np.array([dcrit[i] for i in indices])
        sorted_conn = np.sort(conn)

        ax2.plot(sorted_conn, sorted_dcrit, '+-', linestyle='dotted',
                 label=r'Nestedness $\approx$' + str(nest), markersize=markersize, markeredgewidth=markeredgewidth)
    ax2.set_title(r'$N_R=$' + str(NR) + r', $N_S=$' + str(NS))
    ax2.set_xlabel(r'Connectance', fontsize=fs)
    ax2.set_ylabel(r'$\Delta^*$', fontsize=fs)
    ax2.legend()

    plt.show()
    plt.clf()
