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
title = r'Nr25_Nc25_Nest0.1_Conn0.0832 (1000 runs per point)'
fs = 12
markeredgewidth = 2
markersize = 6
save_string = save_folder + '/' + save_name + '.pdf'
save_string = save_folder + '/proba_geq_one_large.pdf'
delta_crit = 0.0250847
scale = 227.485


def sigmoidal(x, a, b):
    result = []
    for el in x:
        result.append(1 / (1 + np.exp(-b * (el - a))) - 0.5)
    return np.array(result)


def polynomial(x, a):
    result = []
    for el in x:
        sum = 0.
        for i in range(len(a)):
            sum += a[i] * (el**i)
        result.append(sum)
    return np.array(result)


# data = np.loadtxt(folder + '/' + filename + '.out')
# delta = data[:, 0]
# proba = data[:, 1]

delta = [0.0157463,0.0174171,0.0190878,0.0207586,0.0224293,0.0241001,0.0257708,0.0274416,0.0291123,0.0307831]
proba = [-0.44,-0.38,-0.18,-0.26,-0.22,1.11022e-16,0.02,0.22,0.1,0.32]
interval = np.linspace(min(delta), max(delta), 1000)
a = [-0.00676963, 12.3858, -4318.25, 567999, -
     2.30935e+07, 3.91162e+08, -2.41507e+09]
fitted_curve = polynomial(interval, a)
fitted_curve = sigmoidal(interval, delta_crit, scale)
# for i in range(len(delta)):
#     print(str(delta[i]) + ',')

fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(delta, proba, '+', linestyle='none', markersize=markersize)
ax1.plot(interval, fitted_curve)
ax1.set_title(title)
ax1.set_xlabel(r"$\Delta$", fontsize=fs)
ax1.set_ylabel(r"Prob(#extinctions $\geq$ 1)", fontsize=fs)
plt.savefig(save_string)
plt.show()
plt.clf()
