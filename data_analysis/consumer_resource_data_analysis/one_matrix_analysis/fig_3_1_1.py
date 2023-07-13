import pandas as pd
import consumer_resource_data_analysis.consumer_resource_data_analysis as cf
import consumer_resource_data_analysis.data_plotting as dp
import matplotlib.pyplot as plt
import numpy as np

filename = '/Users/Shared/Master/Master_Thesis/results/one_matrix_data_NR25_NS25_verbose-level=1.out'
savename = 'fig_3_1_1'

columns = ['G-matrix path', 'A-matrix path','connG', 'nestG','connA','nestA',
        'alpha_mode','alpha0', 'gamma0','S0','Nsimuls','feasibility probability',
        'stability probability','instability probability','marginality probability',
        'average dominant eigenvalue']

alphamodes = ['fully_connected']
savefolder = '/Users/Shared/Master/Master_Thesis/text/figures/paper'
alpha0_val = 0

### Chargement valeur calculée ####
data = pd.read_csv(filename, comment='#',delimiter=" ", names=columns)

### Calcul théorique ###
maxcoldegG = cf.compute_largest_column_degree("../matrices/Nr25_Nc25/consumption/RandTrix_Nr25_Nc25_Nest0.55_Conn0.2304.txt")
gamma0s = np.linspace(1/maxcoldegG, 1, 100)
th_S0s = cf.theoretical_S0_feasible_boundary(maxcoldegG, gamma0s)

### Figure de comparaison ###
fig, axs = dp.plot_parameters_fixed_alpha0(data, 'feasibility probability', alpha0_val, alphamodes)
axs[0].plot(gamma0s, th_S0s, linestyle='dashed', color='black', marker=None, markersize=0, label="Theoretical threshold")
axs[0].set_title('Feasibility probability')
fig.tight_layout()
fig.savefig(savefolder+'/'+savename+'.pdf')
plt.show()
