import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy.optimize import curve_fit

save_folder='plots/Mon10Feb'
data_folder = 'data_output'
filename = 'test'
title = r'$\alpha_0 = 0, \ \sigma_0 = 0.5, S_0=0.1$ (NR=NS=25 fully connected)'

def power_fit(x, a, b):
    return b*np.power(x-1,a)

data=np.loadtxt(data_folder+'/'+filename+'.out')
S0 = data[:,0]
gamma0 = data[:,1]
unstable_full = data[:, 2]
marginal_full = data[:, 3]
stable_full = data[:, 4]
unstable_eff = data[:, 5]
marginal_eff = data[:, 6]
stable_eff = data[:, 7]

fig1=plt.figure(1)
ax1=fig1.add_subplot(311)
ax1.set_title(title)
ax1.plot(gamma0, unstable_full, marker='', linestyle='solid', label='full')
ax1.plot(gamma0, unstable_eff, marker='',linestyle='solid', label='eff')
ax1.legend()
ax1.set_ylabel('Unstable')

ax2=fig1.add_subplot(312)
ax2.plot(gamma0, marginal_full, marker='', linestyle='solid')
ax2.plot(gamma0, marginal_eff, marker='',linestyle='solid')
ax2.set_ylabel('Marginal')


ax3=fig1.add_subplot(313)
ax3.plot(gamma0, stable_full, marker='', linestyle='solid')
ax3.plot(gamma0, stable_eff, marker='',linestyle='solid')
ax3.set_xlabel(r'$\gamma_0$')
ax3.set_ylabel('Stable')

fig1.tight_layout()
fig1.savefig(save_folder+'/'+filename+'.pdf')

plt.show()
