import glob
import numpy as np

matrix_folder_path = "optimal_matrices/syntrophy/optimal_LRI"

for a in glob.glob(matrix_folder_path + "/RandTrix*.txt"):
    print(a)
