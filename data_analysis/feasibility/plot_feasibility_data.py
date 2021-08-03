import consumer_resource_data_analysis as cf
import matplotlib.pyplot as plt

feasible_volume=dict({
    'intrashift': 0.1,
    'intershift': 0.5,
    'width': 1,
    'data file': 'data_output/all_mat_feasible_volume.csv',
    'save name': 'main_figures/feasibility/feasible_volume_Nr25_Nc25.pdf',
    'alpha mode': ['optimal_matrix','fully_connected', 'random_structure'],
    'alpha0': cf.alpha0
})

feasible_decay_rates = dict({
    'data file': 'data_output/all_mat_feasible_decay_rates.csv',
    'save name': 'main_figures/feasibility/feasible_decay_rates',
    'alpha mode': ['optimal_matrix','fully_connected', 'random_structure'],
})


figures_to_plot = ['feasible volume', 'feasible decay rate']


cf.plot_data(figures_to_plot, feasible_volume, feasible_decay_rates, type='feasible')
