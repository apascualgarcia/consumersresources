import consumer_resource_data_analysis as cf
import matplotlib.pyplot as plt

data_file="data_output/all_mat_local_dynamical_stability.csv"
save_file='plots/lds_volume_Nr25_Nc25.pdf'
cf.alpha_mode=['optimal_matrix', 'fully_connected', 'random_structure']


fig = plt.figure(1)
ax = fig.add_subplot(111)
intrashift= 0.1
intershift= 0.5


shift = [intershift, intrashift]
width=1

ax = cf.plot_lds_volume(ax, data_file, width, shift, cf.alpha0, cf.alpha_mode)

fig.tight_layout()
fig.savefig(save_file)
