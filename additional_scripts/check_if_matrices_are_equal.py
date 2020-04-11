import numpy as np
import os

folder1='optimal_matrices/consumption/Nr25_Nc25'
folder2='matrices/Nr25_Nc25'

files1 = [f for f in os.listdir(folder1) if not f.startswith('.') and f.endswith('.txt') and f.startswith('RandTrix')]
files2 = [f for f in os.listdir(folder2) if not f.startswith('.') and f.endswith('.txt')and f.startswith('RandTrix')]

for i in range(len(files1)):
    mat1 = np.loadtxt(folder1+'/'+files1[i])
    mat2 = np.loadtxt(folder2+'/'+files2[i])

    if not (mat1==mat2).all():
        print('Careful, mat '+ files1[i]+' is not equal to its optimal form')
        print(mat1-mat2)
