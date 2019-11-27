import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
# rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
# # for Palatino and other serif fonts use:
# # rc('font',**{'family':'serif','serif':['Palatino']})
# rc('text', usetex=False)

folder = './data_output'
filename = 'test'
save_folder = './plots'
save_name = filename
title = r'Nr25_Nc25_Nest0.3_Conn0.1296 (100 runs per point)'
fs = 12
markeredgewidth = 2
markersize = 6
save_string = save_folder + '/' + save_name + '.pdf'
save_string = save_folder + '/proba_geq_one_large.pdf'

data = np.loadtxt(folder + '/' + filename + '.out')
delta = data[:, 0]
proba = data[:, 1]


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(delta, proba, '+', linestyle='none', markersize=markersize)
ax1.set_title(title)
ax1.set_xlabel(r"$\Delta$", fontsize=fs)
ax1.set_ylabel(r"Prob(#extinctions $\geq$ 1)", fontsize=fs)
plt.savefig(save_string)
plt.clf()
