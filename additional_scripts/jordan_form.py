import numpy as np
from sympy import Matrix

filename='matrices/Nr25_Nc25/RandTrix_Nr25_Nc25_Nest0.4_Conn0.128.txt'
gamma = np.loadtxt(filename)
print(np.linalg.eigvals(gamma))

m = Matrix(gamma)
P, J = m.jordan_form()
synt = np.transpose(gamma)
for i in range(len(gamma)):
    for j in range(len(gamma[i])):
        if(gamma[i,j]>0):
            synt[j,i]=0
print(synt)
