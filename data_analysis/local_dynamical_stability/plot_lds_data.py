import consumer_resource_data_analysis as cf
import matplotlib.pyplot as plt

ldstable_volume=dict({
    'intrashift': 0.1,
    'intershift': 0.5,
    'width': 1,
    'data file': 'data_output/all_mat_dynamically_stable_volume.csv',
    'save name': 'plots/lds_volume_Nr25_Nc25.pdf',
    'alpha mode': ['optimal_matrix', 'fully_connected', 'random_structure'],
    'alpha0': cf.alpha0
})

ldstable_decay_rates = dict({
    'file name': 'data_output/all_mat_dynamically_stable_decay_rates.csv',
    'save name': 'plots/lds_decay_rates',
    'alpha mode': ['optimal_matrix', 'fully_connected', 'random_structure'],
})


figures_to_plot = ['ld stable volume', 'ld stable decay rate']


cf.plot_data(figures_to_plot, ldstable_volume, ldstable_decay_rates, type='ld stable')
