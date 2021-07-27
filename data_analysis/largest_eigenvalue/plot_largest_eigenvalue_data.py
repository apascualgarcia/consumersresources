import consumer_resource_data_analysis as cf
import matplotlib.pyplot as plt

data_file="data_output/all_mat_dominant_eigenvalue.csv"
save_file='plots/dominant_eigenvalue_Nr25_Nc25.pdf'
cf.alpha_mode=['fully_connected', 'random_structure']


fig = plt.figure(1)
ax = fig.add_subplot(111)
intrashift= 0.1
intershift= 0.5


shift = [intershift, intrashift]
width=1


ax = cf.plot_largest_eigenvalue(ax, data_file, width, shift, cf.alpha0, cf.alpha_mode)
ax.set_yscale('log')
fig.tight_layout()
fig.savefig(save_file)
