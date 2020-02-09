import fnmatch
from os import listdir
from os.path import isfile, join

save_folder = 'data_output/'
filename='critical_delta_convergence_study_conv=1e-1_configuration_comparison_NR25_NS25_s05_a0_matrix_list_NR25_NS25_test_S0=0.1_gamma0=0.5'
potential_candidates = [f for f in listdir(save_folder) if isfile(join(save_folder, f))]

# CHANGE POTENTIAL WILDCARD HERE
file_list=fnmatch.filter(potential_candidates, '*.out')

total_strings = []
for file in file_list:
    with open(save_folder+'/' + file ) as f:
        toapp = f.read().splitlines()
        for str in toapp:
            total_strings.append(str)

filename_comb = filename + '_combined'
f = open(save_folder+'/' + filename_comb + '.out', 'w+')
for str in total_strings:
    f.write(str + '\n')
print("Combined files")
f.close()
