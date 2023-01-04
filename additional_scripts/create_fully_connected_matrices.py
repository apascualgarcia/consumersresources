import numpy as np
import os
N = range(1,101)
for n in N:
    folder = 'matrices/Nr'+str(n)+'_Nc'+str(n)
    if not os.path.exists(folder):
        os.makedirs(folder)
    matrix = np.ones((n,n), dtype=int)
    np.savetxt(folder+'/RandTrix_Nr'+str(n)+'_Nc'+str(n)+'_Nest1_Conn1.txt', matrix, fmt='%d')
