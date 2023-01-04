import numpy as np
from os import listdir
from os.path import isfile, join
import re

matrix_folder="test_matrices"
matrices=[matrix_folder+'/'+f for f in listdir(matrix_folder) if isfile(join(matrix_folder, f)) and not f.startswith('.')]
for mat in matrices:
    [nest, conn] = re.findall("\d+.\d+", mat)
    new_name = mat[:mat.find(nest)]+str(round(float(nest),2)).rstrip('0')+mat[mat.find(nest)+len(nest):-4].strip('0')+mat[-4:]
    np.savetxt(new_name, np.loadtxt(mat), fmt='%d')
