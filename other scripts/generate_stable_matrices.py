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
import networkx as nx
#from networkx.algorithms import bipartite

def build_gamma_matrix(NS_, NR_, param):
    gamma = np.zeros((NS_, NR_))
    if param=='consumer specialist' or param=='mainly specialist':
        eps = 1e-1
        for i in range(min(NS_, NR_)):
            gamma[i,i] = np.random.normal(1., 0.05)
        if param == 'mainly specialist':
            for i in range(NS_):
                for j in range(NR_):
                    if i!=j and np.random.choice(2)==0:
                        gamma[i,j] = np.random.uniform(low=0., high=0.5)
    if param=='nested system' or param=='mainly specialist nested':
        for i in range(min(NS_, NR_)):
            gamma[i,i] = np.random.normal(1., 0.05) 
            if param=='nested system':
                maxvalue = gamma[i,i]
                for j in range(len(gamma[i])):
                    if i < j:
                        gamma[i,j]=np.random.uniform(low=0., high=maxvalue)
                        maxvalue=gamma[i,j]
                        
            elif param=='mainly specialist nested':
                eps = 1e-1
                for j in range(len(gamma[i])):
                    if i < j:
                        gamma[i,j]=np.random.uniform(low=0., high=eps)                
    
    
    return gamma
            
def build_biological_Gamma(NR_, NS_):
    # draw the column sum
    eps = np.random.uniform(low=-0.5, high=0.0, size=NS_)
    
    mincolumnvalue = -0.5
       
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
    
    
def build_beta_Gamma(NR_, NS_, param):
    eps = 0
    beta = np.zeros((NS_, NR_))
    
    for i in range(len(beta)):
        for j in range(len(beta[i])):
            if i==j:
                beta[i,i] = np.random.normal(1., 0.05)
                
    Gamma = np.zeros((NR_, NS_))
    for i in range(len(Gamma)):
        for j in range(len(Gamma[i])):
            if i==j:
                Gamma[i,i] = np.random.normal(-1., 0.05)
    return (beta, Gamma)
    
def diagdom(a, param):
    (minD, maxD) = param
    for i in range(len(a)):
        if minD*abs(a[i,i]) < maxD*sum([abs(a[i,j]) for j in range(len(a[i])) if j!=i]):
            return False
    return True
def prod_strong_dd(D, NR_, NS_):
    minD = np.min(D)
    maxD = np.max(D)
    (beta, Gamma) = build_beta_Gamma(NR_, NS_, (minD, maxD))
    prod = np.dot(beta, Gamma)
    while(not(diagdom(0.5*(prod+np.transpose(prod)), (minD, maxD)))):
        (beta, Gamma) = build_beta_Gamma(NR_, NS_, (minD, maxD))
        prod = np.dot(beta, Gamma)
    return (beta, Gamma)
def draw_parameters(NR_, NS_, param):
    (mode, gamma_mode, degree_of_dissipation)=param
    # the idea is to draw beta, Gamma and D such that the stability condition will be easily satisfied
    # and from there reconstruct the rest of the model
    Seq = np.random.uniform(low=0, high=1., size=NS_)
    Req = np.random.uniform(low=0., high=1., size=NR_)
    # take uniform sigma
    sigma = np.random.uniform(low=0.5, high=0.5, size=(NS_, NR_))
    D = np.random.uniform(low=0., high=1., size=NR_)
    
    if mode=='strong dd':
        # generate two matrices such that their product is strongly diagonally dominant
        (beta, Gamma) = prod_strong_dd(D, NR_, NS_)
        # from the parameters drawn before, we rebuild the entire system
        gamma = np.zeros(shape=(NS_, NR_))
        for i in range(NS_):
            for j in range(NR_):
                gamma[i,j]=beta[i,j]/(Seq[i]*sigma[i,j])

        mu = D-np.dot(np.transpose(gamma), Seq)
        alpha = Gamma-np.dot(np.diag(Req), np.transpose(gamma))
        
    elif mode == 'biological system':
        beta = np.random.uniform(low=0., high=1., size=(NS_, NR_))
        Gamma = build_biological_Gamma(NR_, NS_)
        
        print("Gamma found, ", Gamma)
        
        # from the parameters drawn before, we rebuild the entire system
        gamma = np.zeros(shape=(NS_, NR_))
        for i in range(NS_):
            for j in range(NR_):
                gamma[i,j]=beta[i,j]/(Seq[i]*sigma[i,j])

        mu = D-np.dot(np.transpose(gamma), Seq)
        alpha = Gamma-np.dot(np.diag(Req), np.transpose(gamma))
    
    elif mode=='gamma_construction_restrictive':
        eps = 0.5
        mu = np.random.uniform(low=0., high=1., size=NR_)
        gamma = build_gamma_matrix(NS_, NR_)
        alpha = np.zeros((NR_, NS_))
        for i in range(NR_):
            for j in range(NS_):
                    alpha[i,j] = np.random.uniform(low=0., high=Req[i]*gamma[j,i])
        Gamma = alpha-np.dot(np.diag(Req), np.transpose(gamma))
        
    elif mode == 'gamma_construction':
        repeat = True
        while repeat:
            eps = 0.5
            mu = np.random.uniform(low=0., high=1., size=NR_)
            gamma = build_gamma_matrix(NS_, NR_, gamma_mode)
            alpha = np.random.uniform(low=.0, high=1., size=(NR_, NS_))
            for i in range(NR_):
                for j in range(NS_):
                    if gamma[j,i]>0:
                        alpha[i,j]=0
            Gamma = alpha-np.dot(np.diag(Req), np.transpose(gamma))
            repeat = not(natural_condition_Gamma(Gamma, degree_of_dissipation))
            if repeat :
                print("system rejected because Gamma does not fulfill energy dissipation")
            
        beta = np.dot(np.diag(Seq), np.multiply(sigma, gamma))
        print('beta=')
        print(beta)
        print(np.linalg.eigvals(beta))
        print('Gamma=')
        print(Gamma)
        print(np.linalg.eigvals(Gamma))
        print('beta*Gamma=')
        print(np.dot(beta, Gamma))
        print(np.linalg.eigvals(np.dot(beta,Gamma)))
        # print('beta*Gamma-Gamma*beta=')
#         print(np.dot(beta, Gamma)-np.dot(Gamma, beta))
    else :
        print('draw_parameters does not know how to implement parameters this way!')

    # find the remaining parameters
    delta = np.dot(np.multiply(sigma, gamma), Req)
    lambdaa = np.dot(np.diag(mu), Req)-np.dot(Gamma, Seq)

    return (Req, Seq, alpha, gamma, sigma, lambdaa, delta, mu)
def natural_condition_Gamma(Gamma, deg):
    if deg=='random':
        for j in range(len(Gamma[0])):
            if sum([Gamma[i,j] for i in range(len(Gamma))]) > 0:
                return False
        return True
    else:
        for j in range(len(Gamma[0])):
            somme = sum([Gamma[i,j] for i in range(len(Gamma))])
            if somme > deg or somme < 2*deg :
                return False
        return True
            

# main function of the code, tries to generate a physical system following the stability criterion
def generate_physical_parameters(NR_, NS_, mode):
    print('Looking for physical system...')
    system = draw_parameters(NR_, NS_, mode)
    repeat = not(physical(system))
    if repeat :
        print("System rejected because it was not physical")
    i=1
    while repeat:
        system = draw_parameters(NR_, NS_, mode)
        i+=1
        repeat = not(physical(system))
        if repeat :
            print("System rejected because it was not physical")
        
    print('Physical system found at iteration ', i)
    return system
def physical(A):
    #print('In debug mode, every system is labelled physical')
    #return True
    for a in A:
        for i in range(len(np.ravel(a))):
            if np.ravel(a)[i] < 0:
                return False
    return True
def generate_jacobian_from_parameters(param):
    (Req, Seq, alpha, gamma, sigma, lambdaa, delta, mu) = param
    NR_ = len(alpha)
    NS_ = len(gamma)
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
    
    

# creates a bipartite food network using the consumption rate matrix
def create_syntrophic_network(param):
    
    (Req, Seq, alpha, gamma, sigma, lambdaa, delta, mu)=param
    
    # the food network
    Food = nx.DiGraph()

    # number of species
    NS_ = len(gamma)
    # number of resources
    NR_ = len(gamma[0])

    resources_nodes = ['R'+str(i+1) for i in range(NR_)]
    species_nodes = ['S'+str(i+1) for i in range(NS_)]
    
    for i in range(len(resources_nodes)):
        Food.add_node(resources_nodes[i], bipartite=0, eq_value=Req[i])
    for i in range(len(species_nodes)):
        Food.add_node(species_nodes[i], bipartite=1, eq_value=Seq[i])
    

    # add edges related to food consumption
    for i in range(NS_):
        for j in range(NR_):
            if gamma[i, j] > 0.:
                Food.add_edge(resources_nodes[j],species_nodes[i], weight=sigma[i,j]*gamma[i,j]*Req[j], edge_type='consumption')
    # add edges related to syntrophy
    for i in range(NR_):
        for j in range(NS_):
            if alpha[i,j]>0.:
                Food.add_edge(species_nodes[j], resources_nodes[i], weight=alpha[i,j]*Seq[j], edge_type='byproduct release')

    return Food
# returns a figure containing the plotting for the desired network
def draw_syntrophic_network(axis, A, network_type = 'consumption'):

    colors = dict({'consumption': 'red', 'byproduct release': 'blue'})
        
    resources_nodes = [u for (u,d) in A.nodes(data=True) if d['bipartite']==0]
    species_nodes = [u for (u,d) in A.nodes(data=True) if d['bipartite']==1]
    
    resources_sizes = [1000*d['eq_value'] for (u,d) in A.nodes(data=True) if d['bipartite']==0]
    species_sizes = [1000*d['eq_value'] for (u,d) in A.nodes(data=True) if d['bipartite']==1]
    
    relevant_edges = [(u,v) for (u,v,d) in A.edges(data=True) if d['edge_type']==network_type]
    relevant_weights = [3*d['weight'] for (u,v,d) in A.edges(data=True) if d['edge_type']==network_type]
    
    pos = nx.drawing.layout.bipartite_layout(A, species_nodes, align='horizontal')
    
    
    
    plt.sca(axis)
    #nx.draw_networkx_labels(A, pos=pos)
    nx.draw_networkx_edges(A, pos=pos, edgelist=relevant_edges, width = relevant_weights, edge_color=colors[network_type])
    nx.draw_networkx_nodes(A, nodelist = resources_nodes, pos=pos, node_size=resources_sizes, node_shape='s')
    nx.draw_networkx_nodes(A, nodelist = species_nodes, pos=pos, node_size=species_sizes, node_shape='h')
    axis.set_title(network_type+' network')
    
    return 
def print_system(system):
    (Req, Seq, alpha, gamma, sigma, lambdaa, delta, mu) = system
    
    print('The model has the following properties : ')
    print('Seq = ', Seq)
    print('delta = ', delta)
    print('Req = ', Req)
    print('mu = ', mu)
    print('lambda = ', lambdaa)  
    print('alpha=')
    print(alpha)
    print('gamma=')
    print(gamma)
    print('sigma=')
    print(sigma)
    
    J, trash = generate_jacobian_from_parameters(system)
    print("The Jacobian at equilibrium is given by : ")
    print(J)
    eigvals = np.sort(np.linalg.eigvals(J))
    print('Its eigenvalues are ', eigvals)
    stable = np.max(np.real(eigvals))<=1e-15
    if stable:
        stable = 'stable'
    else:
        stable = 'unstable'
    print('System is', stable)
    
    return 
    
    
    


# def generate_appropriate_matrices(NS_, NR_):
#     gamma = np.zeros((NS_, NR_))
#     alpha = np.zeros((NR_, NS_))
#     for i in range(min(NS_, NR_)):
#         gamma[i,i] = np.random.normal(1., 0.05)
#         alpha[i,i] = np.random.normal(1., 0.05)
#     return gamma, alpha
# def stability_criterion(beta, Gamma, D):
#     prod = np.dot(beta, Gamma)
#     prod = 0.5*(prod+np.transpose(prod))
#
#     minimum = np.min(D)
#     maximum = np.max(D)
#
#     print(prod)
#
#     for i in range(len(prod)):
#         somme = 0.
#         for j in range(len(prod[i])):
#             if i !=j:
#                 somme+=abs(prod[i,j])
#         if minimum*abs(prod[i,i])<maximum*somme:
#             return False
#
#     return True
# def create_stable_system(NR_, NS_):
#     (Req, Seq, alpha, gamma, sigma, lambdaa, delta, mu) = generate_stable_parameters(NR_, NS_, 'default')
#     D = mu + np.dot(np.transpose(gamma), Seq)
#     D = np.diag(D)
#     beta = np.dot(np.diag(Seq), np.multiply(sigma, gamma))
#     Gamma = alpha - np.dot(np.diag(Req), np.transpose(gamma))
#
#
#     while(not(stability_criterion(beta, Gamma, D))):
#         (Req, Seq, alpha, gamma, sigma, lambdaa, delta, mu) = generate_physical_parameters(NR_, NS_, 'default')
#         D = mu + np.dot(np.transpose(gamma), Seq)
#         D = np.diag(D)
#         beta = np.dot(np.diag(Seq), np.multiply(sigma, gamma))
#         Gamma = alpha - np.dot(np.diag(Req), np.transpose(gamma))
#
#
#     return (Req, Seq, alpha, gamma, sigma, lambdaa, delta, mu)
    
NR = 3
NS = 3
dissipation = 'random'
mode = ('gamma_construction', 'mainly specialist nested', dissipation)
system = generate_physical_parameters(NR, NS, mode)
print_system(system)
network = create_syntrophic_network(system)

fig, (ax1, ax2) = plt.subplots(1,2)
draw_syntrophic_network(ax1, network, 'consumption')
draw_syntrophic_network(ax2, network, 'byproduct release')
plt.tight_layout()
plt.show()