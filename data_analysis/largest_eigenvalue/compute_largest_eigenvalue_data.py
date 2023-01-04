import consumer_resource_data_analysis as cf



filename = 'largest_eigenvalue/opt_mat_largest_eigenvalue_NR25_NS25_100_points_full_rank_opt_consumption_mat_NR25_NS25'
optimal_LRI_folder='9Jul21'
consumption_matrix_folder='optimal_matrices/consumption/Nr25_Nc25'
save_file = 'data_output/opt_mat_dominant_eigenvalue.csv'
to_compute = ['av. dominant eigenvalue']
alpha_mode=['optimal_matrix']


cf.compute_largest_eigenvalue_data(alpha_mode, cf.alpha0, filename, optimal_LRI_folder, consumption_matrix_folder, to_compute, save_file)
