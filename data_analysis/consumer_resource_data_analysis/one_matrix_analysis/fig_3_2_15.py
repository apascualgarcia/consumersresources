import pandas as pd
import consumer_resource_data_analysis as cf
import consumer_resource_data_analysis.data_plotting as dp
import matplotlib.pyplot as plt


filename = '/Users/Shared/Master/Master_Thesis/results/one_matrix_data_NR25_NS25_verbose-level=1.out'
savename = 'fig_3_2_15'

columns = ['G-matrix path', 'A-matrix path','connG', 'nestG','connA','nestA',
        'alpha_mode','alpha0', 'gamma0','S0','Nsimuls','feasibility probability',
        'stability probability','instability probability','marginality probability',
        'average dominant eigenvalue']

alphamodes = ['fully_connected', 'random_structure', 'optimized_matrix']
savefolder = '/Users/Shared/Master/Master_Thesis/text/figures/paper'

data = pd.read_csv(filename, comment='#',delimiter=" ", names=columns)

connG = pd.unique(data['connG'])[0]
nestG = pd.unique(data['nestG'])[0]
alpha0s = sorted(pd.unique(data['alpha0']))


extra_params = [connG, nestG]

colors = [dp.cmap(i/(len(dp.alpha0))) for i in range(len(dp.alpha0))]

fig, axs = dp.plot_levels(data, colors, 'stability probability', alphamodes)
fig.savefig(savefolder+'/'+savename+'.pdf')
plt.show()
