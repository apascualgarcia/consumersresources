import consumer_resource_data_analysis as cf


filename = 'feasibility/all_mat_feasibility_NR25_NS25_100_points_full_rank_opt_consumption_mat_NR25_NS25'
optimal_LRI_folder='optimal_LRI_Nr25_Nc25'
consumption_matrix_folder='optimal_matrices/consumption/Nr25_Nc25'
save_file = 'data_output/test.csv'
to_compute = ['feasibility volume']
alpha_mode=['fully_connected', 'random_structure', 'no_release_when_eat']


# filename = 'feasibility/feasibility_NR50_NS25_full_rank_opt_consumption_mat_NR50_NS25'
# optimal_LRI_folder='optimal_LRI_Nr50_Nc25'
# consumption_matrix_folder='optimal_matrices/consumption/Nr50_Nc25'
# matrix_set='S_{50}'
cf.compute_feasibility_data(alpha_mode, cf.alpha0, filename, optimal_LRI_folder, consumption_matrix_folder, to_compute, save_file)
