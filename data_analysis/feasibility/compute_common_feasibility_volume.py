import numpy as np
import sys
from common_functions import remove_strings_from_file
import copy

filename = sys.argv[1]
mat_folder=sys.argv[2]
save_name = sys.argv[3]

remove_strings_from_file(mat_folder, filename)
data = np.loadtxt(filename+'_filtered.out')

NR=data[:,0]
NS=data[:,1]
nestedness=data[:,2]
connectance=data[:,3]
gamma0=data[:, 4::3]
S0=data[:, 5::3]
feasibility=data[:, 6::3]

# COMPUTE COMMON FULL FEASIBILITY VOLUME : VOLUME where all of them are fully feasible
full_feas_indices=[j for j in range(len(feasibility[0])) if feasibility[0,j]==1.]
for i in range(1,len(data)):
    oldffi=copy.deepcopy(full_feas_indices)
    full_feas_indices=[j for j in range(len(feasibility[i])) if (feasibility[i,j]==1. and j in oldffi)]
full_feas_points=[[gamma0[0][j], S0[0][j]] for j in full_feas_indices]


np.savetxt(save_name, full_feas_points)
