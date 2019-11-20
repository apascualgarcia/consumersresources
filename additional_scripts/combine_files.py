filename = 'critical_delta_NR25_NS25'
file_list = [filename + '_' + str(i) for i in range(1, 7)]

total_strings = []
for file in file_list:
    with open('./data_output/' + file + '.out') as f:
        toapp = f.read().splitlines()
        for str in toapp:
            total_strings.append(str)

filename_comb = filename + '_combined'
f = open('./data_output/' + filename_comb + '.out', 'w+')
for str in total_strings:
    f.write(str + '\n')
f.close()
