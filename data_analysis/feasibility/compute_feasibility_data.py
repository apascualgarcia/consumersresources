import consumer_resource_data_analysis as cf

volume_data=dict({
    'file name': 'feasibility/opt_mat_feasibility_NR25_NS25_test_run_full_rank_opt_consumption_mat_NR25_NS25',
    'OM folder': '9Jul21',
    'G matrices folder': 'optimal_matrices/consumption/Nr25_Nc25',
    'save file': 'data_output/opt_mat_feasible_volume.csv',
    'alpha_mode': ['optimal_matrix'],
    'alpha0': cf.alpha0
})

decay_rate_data=dict({
    'file name': 'feasibility/opt_mat_feasibility_NR25_NS25_test_run_full_rank_opt_consumption_mat_NR25_NS25',
    'OM folder': '9Jul21',
    'G matrices folder': 'optimal_matrices/consumption/Nr25_Nc25',
    'save file' : 'data_output/opt_mat_feasible_decay_rates.csv',
    'alpha_mode': ['optimal_matrix'],
    'alpha0_range': cf.alpha0
})

to_compute = ['feasible volume', 'feasible decay rate']

cf.compute_feasibility_data(to_compute, volume_data, decay_rate_data)
