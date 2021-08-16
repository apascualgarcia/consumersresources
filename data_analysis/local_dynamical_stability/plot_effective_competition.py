import consumer_resource_data_analysis as cf
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

columns=['feasible volume', 'stable volume', 'unstable volume', 'marginal volume',
                                'rate of return', 'mean(C)', 'dom C eigval', 'trace(C)', 'dom B eigval']
data_file = 'data_output/eff_comp_all_data_NR25_NS25_full_rank_opt_consumption_mat_NR25_NS25_'
data_suffix = '_9Jul21_verbose-level=1_alpha0=0.out'

alpha_mode = ['optimal_matrix', 'fully_connected', 'random_structure']

x_plot = 'dom B eigval'
y_plot = 'stable volume'
y_scale = 'linear'
save_name = 'domBeigval_stability'

fig = plt.figure()
ax = fig.add_subplot(111)

for a in alpha_mode:
    df = pd.DataFrame(data=np.loadtxt(data_file+a+data_suffix, usecols=[1,2,3,4,5,6,7,8,9]), columns=columns)
    df[y_plot]=df[y_plot].abs()
    df.plot(x=x_plot, y=y_plot, ax=ax, linestyle='', label=cf.alpha_mode_label[a], color=cf.alpha_mode_colours[a])

ax.set_ylabel(y_plot)
ax.legend()
ax.set_yscale(y_scale)
fig.tight_layout()
fig.savefig('plots/'+save_name+'.png', dpi=200)
