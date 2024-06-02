#============ VARIABLES TO GENERATE FIGURE 2 =======================#
plotting_properties=dict({
    'intrashift': 0.1,  # distance between two boxes at the same alpha0
    'intershift': 0.5,  # distance between two boxes (last at alpha0 and first at the next value of alpha0)
    'width': 1,         # box width
    'data file': 'results/all_data_NR25_NS25_new_opt_verbose-level=1.out', # location of the data file produced with compute_feasibility_data.py
    'save folder': 'main_figures',
    'mode': 'plot'      # test mode for testing stuff out, plot mode for plotting data
})

#============ DO NOT MODIFY BELOW ==================================#
import consumer_resource_data_analysis.data_loading as dl
import consumer_resource_data_analysis.data_plotting as dp
import matplotlib.pyplot as plt
import numpy as np


plotting_properties['alpha0']=dp.alpha0
data_frame = dl.load_data_frame(plotting_properties['data file'])

plotting_properties['alpha mode'] = data_frame['A-mode'].unique()
N_alphamodes = len(plotting_properties['alpha mode'])

for result in ['feasible volume', 'dynamically stable volume', 'av. dominant eigenvalue']:
    fig, axs = plt.subplots(N_alphamodes, 2, figsize=(3.5*N_alphamodes, 10))
    dp.plot_against_matrix_properties(axs, data_frame,plotting_properties, result)
    fig.tight_layout()
    fig.savefig(plotting_properties['save folder']+"/"+result+".png")