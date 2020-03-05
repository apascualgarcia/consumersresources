import glob
import numpy as np

matrix_folder_path = "./matrices/Nr25_Nc50"

for a in glob.glob(matrix_folder_path + "/RandTrix*.txt"):
    print(a)
