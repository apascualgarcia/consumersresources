
to_modify_path =  'data_output/max_eigenvalues_variable_syntrophy_full_rank_opt_consumption_mat_NR25_NS25_no_release_when_eat_optimal_LRI.out'
modify_with_path = 'matrix_list/full_rank_opt_consumption_mat_NR25_NS25.in'


# read modifier to add in front
to_add=[]
with open(modify_with_path) as fp:
    line=fp.readline()
    while line:
        to_add.append(line.rstrip())
        line=fp.readline()

# open original file
with open(to_modify_path) as f_in:
    lines = [line.rstrip() for line in f_in] # All lines including the blank ones
    lines = [line for line in lines if line] # Non-blank lines
    file_to_modify=lines
i=0
new_file=[]
# replace requested lines
for line in file_to_modify:
    if(not(line[0]=='#')):
        new_file.append(to_add[i]+' '+line)
        i+=1
    else:
        new_file.append(line)

# rewrite original file
with open(to_modify_path,'w') as f:
    for line in new_file:
        f.write(line+'\n')
