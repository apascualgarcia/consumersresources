import consumer_resource_data_analysis as cf

# Three important metrics can be computed from the dynamical stability data:
# the stability volume, decay rate and p-value
to_compute = ['ld stable volume', 'ld stable decay rate', 'ld stable p-value']

volume_data=dict({
    'file name': 'local_dynamical_stability/opt_mat_local_dynamical_stability_NR25_NS25_test_run_full_rank_opt_consumption_mat_NR25_NS25',
    'OM folder': '9Jul21',
    'G matrices folder': 'optimal_matrices/consumption/Nr25_Nc25',
    'save file': 'data_output/opt_mat_dynamically_stable_volume.csv',
    'alpha_mode': ['optimal_matrix'],
    'alpha0': cf.alpha0
})

decay_rate_data=dict({
    'file name': 'local_dynamical_stability/opt_mat_local_dynamical_stability_NR25_NS25_test_run_full_rank_opt_consumption_mat_NR25_NS25',
    'OM folder': '9Jul21',
    'G matrices folder': 'optimal_matrices/consumption/Nr25_Nc25',
    'save file' : 'data_output/opt_mat_dynamically_stable_decay_rates.csv',
    'alpha_mode': ['optimal_matrix'],
    'alpha0_range': cf.alpha0
})

####### DO NOT MODIFY BELOW THIS LINE ################
cf.compute_lds_data(to_compute, volume_data, decay_rate_data)
