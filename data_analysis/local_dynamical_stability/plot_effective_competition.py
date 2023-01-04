import consumer_resource_data_analysis as cf
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

## Plot effective competition figure

data_file = 'data_output/eff_comp_all_data_NR25_NS25_full_rank_opt_consumption_mat_NR25_NS25_'
data_suffix = '_Metamatrices_verbose-level=1_alpha0=0.out' 
alpha_mode = ['optimal_matrix', 'fully_connected', 'random_structure']
x_plot = 'ratio inter- intraspecific competition'
y_plot = 'stable volume'
y_scale = 'linear'
save_name = 'ratio_inter_intra'


######### DO NOT MODIFY BELOW THIS LINE #########
columns=['feasible volume', 'stable volume', 'unstable volume', 'marginal volume',
                                'rate of return', 'mean(C)', 'dom C eigval', 'trace(C)',
                                 'dom B eigval', 'ratio inter- intraspecific competition']
fig = plt.figure()
ax = fig.add_subplot(111)

for a in alpha_mode:
    df = pd.DataFrame(data=np.loadtxt(data_file+a+data_suffix, usecols=range(1, len(columns)+1)), columns=columns)
    if y_scale=='log':
        df[y_plot]=df[y_plot].abs()
    df.plot(x=x_plot, y=y_plot, ax=ax, linestyle='', label=cf.alpha_mode_label[a], color=cf.alpha_mode_colours[a])
ax.set_ylabel(y_plot)
ax.legend()
ax.set_yscale(y_scale)
fig.tight_layout()
fig.savefig('plots/'+save_name+'.png', dpi=200)
