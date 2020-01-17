import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
# rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
# # for Palatino and other serif fonts use:
# # rc('font',**{'family':'serif','serif':['Palatino']})
# rc('text', usetex=False)

folder = './data_output'
filename = 'average_extinctions_large'
save_folder = './plots'
save_name = filename
title = r'Nr25_Nc25_Nest0.3_Conn0.1296 (1000 runs per point)'
fs = 12
error_bar_width = 1.
cap_width = 1.5
markeredgewidth = 2
markersize = 6
errorbar = 'no_errorbar'

data = np.loadtxt(folder + '/' + filename + '.out')
delta = data[:, 0]
av_ext = data[:, 1]
std_dev = data[:, 2]


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
if(errorbar == 'errorbar'):
    ax1.errorbar(delta, av_ext, yerr=std_dev, fmt='+', linestyle='none',
                 markeredgewidth=markeredgewidth, markersize=markersize,
                 elinewidth=error_bar_width, capsize=cap_width)
elif(errorbar == 'no_errorbar'):
    ax1.plot(delta, av_ext, '+', linestyle='none', markersize=markersize)
ax1.set_title(title)
ax1.set_xlabel(r"$\Delta$", fontsize=fs)
ax1.set_ylabel(r"average extinctions", fontsize=fs)
plt.savefig(save_folder + '/' + save_name + '_' + errorbar + '.pdf')
plt.clf()
