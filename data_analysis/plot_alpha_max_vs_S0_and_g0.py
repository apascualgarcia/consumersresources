import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


data_folder='./data_output'
plot_folder='./plots/Mon10Feb'

filename='feasible_alpha_max_S0=0.1_gamma0=0.1_s05_a0'
title=r'25x25 $\gamma$ matrix with $\kappa=0.13$ and $\eta=0.3$, $S_0=0.1$'

def linear_fit(x, a, b):
    return a*x+b

file = data_folder+'/'+filename

data=np.loadtxt(file+'.out')
S0 = data[:,0]
gamma0 = data[:,1]
max_alpha0 = data[:,2]

valid_indices = ~(np.isnan(S0) | np.isnan(gamma0) | np.isnan(max_alpha0))

S0 = S0[valid_indices]
gamma0 = gamma0[valid_indices]
max_alpha0 = max_alpha0[valid_indices]

S0_target = 0.1
target_indices = [i for i in range(len(S0)) if S0[i]==S0_target]
gamma0_target = gamma0[target_indices]
max_alpha0_target = max_alpha0[target_indices]


popt,trash = curve_fit(linear_fit, gamma0_target, max_alpha0_target)
fitted_curve = np.array([linear_fit(g, *popt) for g in gamma0_target])

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
ax1.plot(gamma0_target, fitted_curve, label='linear fit (slope {:2f}'.format(popt[0])+')', linestyle='-', marker='')
ax1.plot(gamma0_target, max_alpha0_target, label='numerical data', marker='o', markersize=2)
ax1.legend()
ax1.set_xlabel(r'$\gamma_0$')
ax1.set_ylabel(r'Maximum feasible $\alpha_0$')
ax1.set_title(title)
fig1.tight_layout()
fig1.savefig(plot_folder+'/'+filename+'.pdf')
plt.show()
