import pandas as pd
import consumer_resource_data_analysis as cf
import consumer_resource_data_analysis.data_plotting as dp
import matplotlib.pyplot as plt


filename = '/Users/Shared/Master/Master_Thesis/results/one_matrix_data_NR25_NS25_verbose-level=1.out'

columns = ['G-matrix path', 'A-matrix path','connG', 'nestG','connA','nestA',
        'alpha_mode','alpha0', 'gamma0','S0','Nsimuls','feasibility probability',
        'stability probability','instability probability','marginality probability',
        'average dominant eigenvalue']

data = pd.read_csv(filename, comment='#',delimiter=" ", names=columns)

connG = pd.unique(data['connG'])[0]
nestG = pd.unique(data['nestG'])[0]
alpha0s = sorted(pd.unique(data['alpha0']))

extra_params = [connG, nestG]

colors=[plt.cm.get_cmap('jet_r')(i/len(alpha0s)) for i in range(len(alpha0s))]

print(sorted(pd.unique(data['gamma0'])))
print(data[data['gamma0']==100])
#dp.plot_levels(data, colors, 'feasibility probability')
#plt.show()
