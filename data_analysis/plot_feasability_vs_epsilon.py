import numpy as np
import matplotlib.pyplot as plt

data_folder='./data_output'
plot_folder='./plots/Mon10Feb'

filename='feasability_vs_epsilon_random_matrix'


file = data_folder+'/'+filename

data=np.loadtxt(file+'.out')
epsilon = data[:,0]
max_alpha0 = data[:,1]

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
ax1.plot(epsilon, max_alpha0)
ax1.set_xlabel(r'$\epsilon$ (relative variance)')
ax1.set_ylabel(r'Maximum feasible $\alpha_0$')
ax1.set_title('Random (low connectance and nestedness) 25x25 matrix')
fig1.tight_layout()
fig1.savefig(plot_folder+'/'+filename+'.pdf')
