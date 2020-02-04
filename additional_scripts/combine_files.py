cores = 4
filename = 'cd_alpha_max_l_0=2_R0=1_configuration_comparison_NR25_NS25_s05_a0_critical_matrix_list_NR25_NS25'
file_list = [filename + '_' + str(i) for i in range(0, cores)]
save_folder = 'data_output/'

total_strings = []
for file in file_list:
    with open(save_folder+'/' + file + '.out') as f:
        toapp = f.read().splitlines()
        for str in toapp:
            total_strings.append(str)

filename_comb = filename + '_combined'
f = open(save_folder+'/' + filename_comb + '.out', 'w+')
for str in total_strings:
    f.write(str + '\n')
print("Combined files")
f.close()
