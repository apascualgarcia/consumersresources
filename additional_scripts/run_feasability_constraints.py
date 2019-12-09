import os
import numpy as np

cmd = 'build/study_feasability_constraints'
config_file = 'config/configuration_comparison_NR25_NS25_s05_a0.in'
alpha_range = np.linspace(0., 1.5, 40)
path_to_save = 'data_output/study_feasability_constraints.out'
to_run = ''
seed = int(0)
for a in alpha_range:
    to_run += cmd + " " + config_file
    to_run += " path_to_save_file="+path_to_save
    to_run += " alpha_0="+ '{0:.5f}'.format(a)
    to_run += " seed_number=" +str(seed)
    to_run += " && "
    seed += 1
to_run = to_run[:-4]
#print(to_run)
os.system(to_run)
