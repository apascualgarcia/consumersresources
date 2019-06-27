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
from mpl_toolkits.mplot3d import Axes3D 
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap


def generate_physical_parameters(NR,NS, type='stability'):
    ## choose physical model
    Req = np.random.uniform(low=0., high=1., size=NR)
    Seq = np.random.uniform(low=0., high=1., size=NS)
    
    
    alpha = np.random.uniform(low=0., high=1., size=(NR,NS))
    gamma = np.random.uniform(low=0., high=1., size=(NS,NR))
    
    if type=='stability':
        for i in range(min(NR, NS)):
            gamma[i,i]= np.random.normal(1., 0.05)
        
    sigma = np.random.uniform(low=0., high=1., size=(NS,NR))
    lambdaa = np.random.uniform(low=0., high=1., size=NR)
    
    delta = np.zeros(NS)
    mu = np.zeros(NR)
    
    for i in range(NR):
        mu[i] = lambdaa[i]
        sumterm = 0.
        for j in range(NS):
            sumterm+=-gamma[j][i]*Req[i]*Seq[j]+alpha[i][j]*Seq[j]
        mu[i] += sumterm
        mu[i] = 1./Req[i]*mu[i]
    
    for i in range(NS):
        for j in range(NR):
            delta[i] += sigma[i][j]*gamma[i][j]*Req[j]
    
    energy_dissipated = False
    
    while(not(physical_parameters(delta, mu)) or not(energy_dissipated)):
        Req = np.random.uniform(low=0., high=1., size=NR)
        Seq = np.random.uniform(low=0., high=1., size=NS)
    
        alpha = np.random.uniform(low=0., high=1., size=(NR,NS))
        gamma = np.random.uniform(low=0., high=1., size=(NS,NR))
        sigma = np.random.uniform(low=0., high=1., size=(NS,NR))
        lambdaa = np.random.uniform(low=0., high=1., size=NR)
        
        if type=='stability':
            for i in range(min(NR, NS)):
                gamma[i,i]= np.random.normal(1., 0.05)
    
        delta = np.zeros(NS)
        mu = np.zeros(NR)
    
        for i in range(NR):
            mu[i] = lambdaa[i]
            sumterm = 0.
            for j in range(NS_):
                sumterm+=-gamma[j][i]*Req[i]*Seq[j]+alpha[i][j]*Seq[j]
            mu[i] += sumterm
    
        for i in range(NS):
            for j in range(NR):
                delta[i] += sigma[i][j]*gamma[i][j]*Req[j]
                
        energy_dissipated = dissipation(alpha, gamma, Req)
        
    return (Req, Seq, alpha, gamma, sigma, lambdaa, delta, mu)
    
def generate_jacobian_from_parameters(param):
    (Req, Seq, alpha, gamma, sigma, lambdaa, delta, mu) = param
    D = np.zeros((NR_,NR_))
    for i in range(NR_):
        D[i][i] = mu[i]
        for j in range(NS_):
            D[i][i] += gamma[j][i]*Seq[j]
            
    G = np.zeros((NR_, NS_))
    for i in range(NR_):
        for j in range(NS_):
            G[i][j]=alpha[i][j]-Req[i]*gamma[j][i]
    
    B = np.zeros((NS_,NR_))
    for i in range(NS_):
        for j in range(NR_):
            B[i][j] = Seq[i]*sigma[i][j]*gamma[i][j]
    
    Z = np.zeros((NS_, NS_))
    for i in range(NS_):
        sumterm = 0.
        for k in range(NR_):
            sumterm+=sigma[i][k]*gamma[i][k]*Req[k]
        Z[i][i] = sumterm - delta[i]

    return np.block([[-D, G], [B, Z]]), (D, G, B, Z)
    

def generate_feasible_jacobian(NR_, NS_):
    return generate_jacobian_from_parameters(generate_physical_parameters(NR_,NS_))
    

def physical_parameter(delta_):
    for d in np.ravel(delta_):
        if d < 0 : 
            return False
    return True

def dissipation(alpha, gamma, Req):
    for i in range(len(alpha)):
        for j in range(len(alpha[0])):
            if alpha[i][j] - gamma[j][i]*Req[i] >= 0:
                return False
    return True

# we want to know when exactly an eigenvalue of J is also an eigenvalue of D    
def find_eigenvalues(J, blocks):
    Jeigvals = np.sort(np.linalg.eigvals(J))
    (D,G,B,Z) = blocks
    Deigvals = np.sort(np.linalg.eigvals(D))
    
    colors = ['blue' for i in range(len(Jeigvals))]
    
    # if eigenvalue is eigenvalue of D, the value is flagged
    for i in range(len(Jeigvals)) :
        for j in range(len(Deigvals)): 
            if abs(Jeigvals[i]-Deigvals[j]) < 1e-3:
                   colors[i]='red'
                   print('Found a D eigenvalue')
               
    return Jeigvals, colors
    
def generate_gaussian_parameters(NR_, NS_, param):
    (gamma0, std, epsilon, heterogeneous_efficiency, eff_value) = param
    
    if heterogeneous_efficiency == True:
        sigma = np.random.uniform(low=0., high=1., size=(NS_,NR_))
    else:
        sigma = np.ones((NS_, NR_))*eff_value
    lambdaa = np.random.uniform(low=0., high=1., size=NR_)
    Req = np.random.uniform(low=0., high=1., size=NR_)
    Seq = np.random.uniform(low=0., high=1., size=NS_)
    
    
    gamma = np.random.normal(loc=gamma0, scale=std, size=(NS_,NR_))
    alpha = np.zeros((NR_, NS_))
    for i in range(NR_):
        for j in range(NS_):
            alpha[i][j]=epsilon*gamma[j][i]*Req[i]
    
    delta = np.zeros(NS_)
    mu = np.zeros(NR_)
    
    
    for i in range(NR_):
        mu[i] = lambdaa[i]
        sumterm = 0.
        for j in range(NS_):
            sumterm+=-gamma[j][i]*Req[i]*Seq[j]+alpha[i][j]*Seq[j]
        mu[i] += sumterm
        mu[i] = 1./Req[i]*mu[i]
    
    for i in range(NS_):
        for j in range(NR_):
            delta[i] += sigma[i][j]*gamma[i][j]*Req[j]
            
    return (Req, Seq, gamma, alpha, sigma, lambdaa, mu, delta)

# generates a jacobian with a gaussianly picked gamma and alpha prop to it
def generate_gaussian_jacobian(NR_, NS_, param):
    
    (Req, Seq, gamma, alpha, sigma, lambdaa, mu, delta) = generate_gaussian_parameters(NR_, NS_, param)
    
    ## choose physical model
    while(not(physical_parameter(delta)) or not(physical_parameter(mu)) or not(physical_parameter(gamma))):
        (Req, Seq, gamma, alpha, sigma, lambdaa, mu, delta) = generate_gaussian_parameters(NR_, NS_, param)
                                
    D = np.zeros((NR_,NR_))
    for i in range(NR_):
        D[i][i] = mu[i]
        for j in range(NS_):
            D[i][i] += gamma[j][i]*Seq[j]
            
    G = np.zeros((NR_, NS_))
    for i in range(NR_):
        for j in range(NS_):
            G[i][j]=alpha[i][j]-Req[i]*gamma[j][i]
    
    B = np.zeros((NS_,NR_))
    for i in range(NS_):
        for j in range(NR_):
            B[i][j] = Seq[i]*sigma[i][j]*gamma[i][j]
    
    Z = np.zeros((NS_, NS_))
    for i in range(NS_):
        sumterm = 0.
        for k in range(NR_):
            sumterm+=sigma[i][k]*gamma[i][k]*Req[k]
        Z[i][i] = sumterm - delta[i]

    return np.block([[-D, G], [B, Z]]), (D, G, B, Z)

# NR = 2
# NS = 3
# Nsimul = 1000
#
# Jeq, blocks = generate_feasible_jacobian(NR, NS)
# eigenvalues, colors = find_eigenvalues(Jeq, blocks)
# for i in range(Nsimul-1):
#     Jeq, blocks = generate_feasible_jacobian(NR, NS)
#     eigvals, cols = find_eigenvalues(Jeq, blocks)
#     eigenvalues = np.concatenate((eigenvalues, eigvals), axis=None)
#     colors = np.concatenate((colors, cols), axis=None)
#
#
# print('Max eigenvalue found : ', max(eigenvalues))
# plt.scatter(np.real(eigenvalues), np.imag(eigenvalues), c= colors, s=2)
# plt.plot([0., 0.], [np.min(np.imag(eigenvalues)), np.max(np.imag(eigenvalues))])
# plt.xlabel(r'Re($\lambda$)')
# plt.ylabel(r'Im($\lambda$)')
# ttle = r'Eigenspectrum for feasible equilibria (Nsimul='+str(Nsimul)+')'
# plt.title(ttle)
# plt.savefig('eigenspectrum_feasible_jacobian.pdf')
# plt.show()

# What simulations show : eigenvalues far from being D eigenvalues, largest eigenvalue always positive but small compared to the intensity of the smallest
# also the spectrum is symmetric around the real axis

# np.random.seed()
#
#
# NR = 2
# NS = 3
# Nsimul = 10000
#
# gamma0 = np.linspace(0.5, 2., 10)
# std=0.05
# epsilon=1.0
# heterogeneous_efficiency = False
# eff_value = 0.1
# #savename = 'eigenspectrum_heterogeneous_eff'
# savename = 'eigenspectrum_homogenous_eff'
# extension = 'png'

# max_eigvals=[]
# for g in gamma0:
#     maxima = []
#     unstable=False
#     i=0
#     while i < Nsimul and not(unstable):
#         param = (g, std, epsilon, heterogeneous_efficiency, eff_value)
#         J, blocks = generate_gaussian_jacobian(NR, NS, param)
#         M = np.max(np.real(np.linalg.eigvals(J)))
#         if  M > 1e-10 :
#             unstable = True
#             print('For gamma0 = ', g, ' found a positive eigenvalue ', M, ' at iteration ', i)
#         i+=1
#
#     max_eigvals.append(unstable)
#
# plt.plot(gamma0, max_eigvals)
# plt.show()

def stability(g, e, param):
    (std, heterogeneous_efficiency, eff_value)=param
    unstable=False
    i=0
    params = (g, std, e, heterogeneous_efficiency, eff_value)
    while i < Nsimul and not(unstable):
        J, blocks = generate_gaussian_jacobian(NR, NS, params)
        M = np.max(np.real(np.linalg.eigvals(J)))
        if  M > 1e-10 :
            unstable = True
        i+=1

    if unstable:
        return 1.
    else :
        return 0.


def run():
    np.random.seed()
    NR = 2
    NS = 3
    Nsimul = 1000

    gamma0 = np.linspace(0.5, 2., 10)
    std=0.05
    epsilon=1.0
    heterogeneous_efficiency = False
    eff_value = 0.1
    #savename = 'eigenspectrum_heterogeneous_eff'
    savename = 'eigenspectrum_homogenous_eff'
    extension = 'png'    
        
    gamma0=np.linspace(0.5, 3., 100)
    epsilon = np.linspace(0., 2., 100)
    param = (std, heterogeneous_efficiency, eff_value)

    z = np.zeros((len(gamma0), len(epsilon)))

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    x = gamma0
    y = epsilon
    X, Y = np.meshgrid(x, y)

    for i in range(len(gamma0)):
        for j in range(len(epsilon)):
            z[i,j]=stability(gamma0[i], epsilon[j], param)
        
    surf = ax.plot_surface(X,Y, np.transpose(z))
    plt.xlabel(r'$\gamma_0$')
    plt.ylabel(r'$\epsilon$')
    plt.title(r'1 means unstable system, 0 means stable')
    plt.tight_layout()
    plt.show()

# we build Gamma such that sum i Gij = eps_j < 0
def build_biological_Gamma(NR_, NS_):
    # draw the column sum
    eps = np.random.uniform(low=-1., high=0., size=NS_)
    print("sum of columns : ", eps)
    
    mincolumnvalue = -1.
       
    Gamma = np.zeros((NR_, NS_))
    for j in range(NS_):
        drawn_value=np.zeros(NR_)
        low_limit = np.random.uniform(low=-mincolumnvalue, high=0.)
        high_limit = eps[j]/NR_
        for i in range(NR_-1):
            drawn_value[i] = np.random.uniform(low=low_limit, high=high_limit)
            low_limit = drawn_value[i]
            high_limit = eps[j]-sum([drawn_value[k] for k in range(i)])/(NR_-i+1)
        drawn_value[NR_-1] = eps[j]-sum([drawn_value[k] for k in range(NR_-1)])
        np.random.shuffle(drawn_value)
        Gamma[:,j] = drawn_value

    return Gamma

def run_Gamma():
    NR, NS = 2,3
    print(build_biological_Gamma(NR, NS))


run_Gamma()
# for j in range(len(gamma0)):
#     eigenvals = []
#     for i in range(Nsimul):
#         param = (gamma0[j], std, epsilon, heterogeneous_efficiency, eff_value)
#         J, trash = generate_gaussian_jacobian(NR, NS, param)
#         eigenvals.append(np.linalg.eigvals(J))
#     eigenvals = np.ravel(np.array(eigenvals))
#     colors = cm.rainbow(np.linspace(0, 1, len(gamma0)))
#     plt.scatter(np.real(eigenvals), np.imag(eigenvals), s=0.5, color=colors[j])
#     plt.plot([0., 0.], [np.min(np.imag(eigenvals)), np.max(np.imag(eigenvals))], color='grey')
#     plt.xlabel(r'Re($\lambda$)')
#     plt.ylabel(r'Im($\lambda$)')
#     plt.title(r'Warmer colour means increasing $\gamma_0$')
#     plt.draw()
# plt.xlim(-1.3,0.1)
# plt.savefig(savename+'.'+extension, format=extension, dpi=400)
# plt.xlim(-0.05, 0.1)
# plt.savefig(savename+'_closeup.'+extension, format=extension, dpi=400)
#