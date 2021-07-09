import numpy as np
import pandas as pd
import consumer_resource_data_analysis as cf
import matplotlib.pyplot as plt

data_file_path = "data_output/optimized_matrix_data_alpha0=1_intra_specific_syntrophy=allowed.csv"


data = pd.read_csv(data_file_path)

fig1 = plt.figure("Nestedness")
ax1 = fig1.add_subplot(111)
xlim = [np.min(data['nestG'])-0.1*(np.max(data['nestG'])-np.min(data['nestG'])), 0.1*(np.max(data['nestG'])-np.min(data['nestG']))+np.max(data['nestG'])]
for connG in cf.all_connectance:
    to_plot = data[data['connG']==connG]
    to_plot.plot(x = 'nestG', y = 'nestA', label=r'$\kappa_G='+str(connG)+'$', ax=ax1, xlim = xlim)
ax1.set_xlabel(r'$\eta_G$')
ax1.set_ylabel(r'$\langle \eta_A \rangle $')
fig1.tight_layout()

fig2 = plt.figure("Connectance")
ax2 = fig2.add_subplot(111)
xlim = [np.min(data['connG'])-0.1*(np.max(data['connG'])-np.min(data['connG'])), 0.1*(np.max(data['connG'])-np.min(data['connG']))+np.max(data['connG'])]
for nestG in cf.all_nestedness:
    to_plot = data[data['nestG']==nestG]
    to_plot.plot(x='connG', y='connA', label=r'$\eta_G='+str(nestG)+'$', ax=ax2, xlim=xlim)
ax2.set_xlabel(r'$\kappa_G$')
ax2.set_ylabel(r'$\langle \kappa_A \rangle$')
fig2.tight_layout()

plt.show()
