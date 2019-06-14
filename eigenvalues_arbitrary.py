#!/usr/local/bin/python3
# import the path for tex to work
import os
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'
import numpy as np
from scipy.optimize import fsolve
#np.random.seed(10)
np.set_printoptions(precision=2)
from numpy.polynomial.polynomial import *
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

# the goal is to find the eigenspectrum of a generic matrix J written as
#      -D G
# J =   B Z  with Z =0


NR = 8
NS = 6
Nsimul = 5000

single_resource=False

def generate_jacobian(NR_,NS_, param):
    (minD, maxD, minG, maxG, minB, maxB) = param
    D = np.diag(np.random.uniform(low=minD, high=maxD, size=NR_))
    G = np.random.uniform(low=minG, high=maxG, size=(NR_,NS_))
    B = np.random.uniform(low=minB, high=maxB, size=(NS_,NR_))
    Z = np.zeros((NS_,NS_))
    
    return np.block([[-D, G],[B, Z]]), (D, G, B, Z)

# if enforce is true, we manually set one random element of G to be equal to the argument given
def generate_special_jacobian(NR_, param, enf):
    jac, (D,G,B,Z) = generate_jacobian(NR_, NR_, param)
    (enforced, value) = enf
    for i in range(NR_):
        for j in range(NR_):
            if i < j:
                G[i][j] = 0
            if i != j:
                B[i][j] = 0
    expected_sol = []
    if enf:
        i = np.random.randint(low=0, high=NR_)
        j = np.random.randint(low=0, high=i+1)
        G[i][j]=value
        
    for i in range(NR_):
        l1 = Z[i][i]-D[i][i]+np.sqrt((D[i][i]+Z[i][i])**2 + 4*B[i][i]*G[i][i]+0j)
        l2 = Z[i][i]-D[i][i]-np.sqrt((D[i][i]+Z[i][i])**2 + 4*B[i][i]*G[i][i]+0j)
        expected_sol.append(0.5*l1)
        expected_sol.append(0.5*l2)
        
    return np.block([[-D, G],[B, Z]]), expected_sol

def generate_feasible_jacobian(NR_, NS_, param):
    Req = np.random.uniform(low=0., high=1., size=NR)
    Seq = np.random.uniform(low=0., high=1., size=NS)
    
    

def generate_eigenspectrum_distribution(NR_, NS_, Nsimul, param):
    eigenspectrum = []
    for i in range(Nsimul):
        J, trash = generate_jacobian(NR, NS, param)
        actual_eig = np.linalg.eigvals(J)
        for j in range(len(actual_eig)):
            eigenspectrum.append(actual_eig[j])
    return np.array(eigenspectrum)

def expected_eigenvalues(*args):
    D, G, B, Z = args
        
    d = np.diagonal(D)
    z = np.diagonal(Z)
    
    nr = len(d)
    ns = len(z)
    
    omega = np.multiply(B, np.transpose(G))
    
    H = (1)
    for k in range(nr):
        H = polymul(H, (d[k],1))
    tilde_H = [(1) for j in range(nr)]
    for j in range(nr):
        for k in range(nr):
            if k!=j:
                tilde_H[j] = polymul(tilde_H[j], (d[j],1))
                
    sols = []
    for i in range(ns):
        sols.append(find_roots(z[i], H,tilde_H,omega[i,:]))
    
    print(sols)
    
    
def find_roots(zi, H, tilde_H, omegai):
    polynom = polymul(H, (-zi, 1))
    for j in range(len(tilde_H)):
        polynom = polysub(polynom, polymul((omegai[j]), tilde_H[j]))
    print(polynom)
    return polyroots(polynom)
    
    
def generate_jacobian_single_res(NS_):
    a = np.random.uniform(low=0., high=1., size=1)
    b = np.random.uniform(low=0., high=1., size=(1,NS_))
    c = np.zeros(shape=(NS_,1))
    cn = np.random.uniform(low=0., high=1., size=1)
    c[-1:] = cn[0]*1.
    d = np.random.uniform(low=-1., high=1., size=NS_-1)
    
    C = cn[0]
    D = b[0][-1]
    
    l1 = a/2.*(np.sqrt(1+4*cn[0]*b[0][-1]/(a**2))-1)
    l2 = -a/2.*(np.sqrt(1+4*cn[0]*b[0][-1]/(a**2))+1)
        
    diag_values = d*1.
    to_add = np.array([l1[0], l2[0]])
    
    d = list(d)+[0.]
    D = np.diag(d)
    
    expected = np.concatenate((diag_values, to_add))

    return np.block([[-a, b], [c, D]]), np.sort(expected)
    
    
def find_max_eigval_special(a,b, nsimul, enf):
    if a <= b:
        eigmax = np.zeros(nsimul)
        for k in range(nsimul):
            param = (0., 1., a, b, 0., 1.)
            J, expected_eigvals = generate_special_jacobian(NR, param, enf)
            eigmax[k] = max(expected_eigvals)
        return np.mean(eigmax)
    else : 
        return np.nan

if NR==1 and single_resource==True:
    print('We generate a random jacobian matrix for the single resource case : ')
    J, expected_eigvals = generate_jacobian_single_res(NS)
    print(J)
    print('The predicted eigenvalues are : ', expected_eigvals)
    print('And the actual eigenvalues are : ', np.sort(np.linalg.eigvals(J)))
# else:
#     minD = 0.
#     maxD = 1.
#     minB = 0.
#     maxB = 1.
#     minG = [-3., -1.5, -.8, -.2]
#     maxG = 0.
#     colors = ['black','red', 'orange', 'yellow']
#     for i in range(len(minG)):
#         eigenspectrum = generate_eigenspectrum_distribution(NR, NS, Nsimul, (minD, maxD, minB, maxB, minG[i], maxG))
#         plt.scatter(np.real(eigenspectrum), np.imag(eigenspectrum), s=0.1, c=colors[i])
#     plt.show()
else : 
    Npoints = 100
    minG = np.linspace(-1., 1., num=Npoints)
    maxG = np.linspace(-1., 1., num=Npoints)
    nsimul = 300
    
    maxeigval = np.zeros((len(minG), len(maxG)))
    
    for i in range(len(minG)):
        for j in range(len(maxG)):
            enf = (False, maxG[j])
            maxeigval[i][j] = find_max_eigval_special(minG[i], maxG[j], nsimul,enf)
    
    f = open('mean_maximum_eigval_wrto_Gamma_NR='+str(NR)+'.txt', "a")
    f.write(str(maxeigval))
    f.close()
    
    rescaled_vals = np.transpose(maxeigval)#*1./(np.max(maxeigval))
    maxvalue = np.max(rescaled_vals[~np.isnan(rescaled_vals)])
    current_map = cm.get_cmap('bwr')
    current_map.set_bad(color='grey')
    plt.imshow(rescaled_vals, cmap = current_map, origin='lower', vmin=-maxvalue, vmax=maxvalue, extent=[minG[0], minG[-1], maxG[0], maxG[-1]])
    cbar = plt.colorbar()
    plt.xlabel('Minimum value of $\Gamma_{ij}$')
    plt.ylabel('Maximum value of $\Gamma_{ij}$')
    cbar.set_label('Maximum eigenvalue', rotation=90)
    plt.savefig('mean_maximum_eigval_wrto_Gamma_NR='+str(NR)+'.pdf')
    plt.savefig('mean_maximum_eigval_wrto_Gamma_NR='+str(NR)+'.png')
    
    plt.show()
    
    
    


