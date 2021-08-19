import consumer_resource_data_analysis as cf
import matplotlib.pyplot as plt

# Two figures about feasibility may be computed :
# feasible volume vs. syntrophy and feasible decay rate vs. syntrophy

figures_to_plot = ['feasible volume', 'feasible decay rate']

feasible_volume=dict({
    'intrashift': 0.1,  # distance between two boxes at the same alpha0
    'intershift': 0.5,  # distance between two boxes (last at alpha0 and first at the next value of alpha0)
    'width': 1,         # box width
    'data file': 'data_output/all_mat_feasible_volume.csv', # location of the data file produced with compute_feasibility_data.py
    'save name': 'main_figures/feasibility/feasible_volume_Nr25_Nc25.pdf',
    'alpha mode': ['optimal_matrix','fully_connected', 'random_structure'],
    'alpha0': cf.alpha0 # list of alpha0 that should be considered
})

feasible_decay_rates = dict({
    'data file': 'data_output/all_mat_feasible_decay_rates.csv',
    'save name': 'main_figures/feasibility/feasible_decay_rates',
    'alpha mode': ['optimal_matrix','fully_connected', 'random_structure'],
})

######## DO NOT MODIFY BELOW THIS LINE #########
cf.plot_data(figures_to_plot, feasible_volume, feasible_decay_rates, type='feasible')
