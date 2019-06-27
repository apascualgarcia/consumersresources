#!/usr/local/bin/python3
# import the path to make graphs with TeX, uncomment if it doesn't work (or put your latex folder location)
import os
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'
import numpy as np
from scipy.optimize import fsolve
from numpy.polynomial.polynomial import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import networkx as nx
from networkx.algorithms import bipartite

# console printing options
np.set_printoptions(precision=2)

# produced graphs display preferences
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


class Model:
    def __init__(self, param='None'):
        self.data = param
    
    # internal function to create the food network that will be drawn later
    def __create_food_network(self):   
        (Req, Seq, alpha, gamma, sigma, lambdaa, delta, mu)=self.retrieve_model_parameters()
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
        self.data['food network'] = Food
        return
      
    # gives back the model parameters
    def retrieve_model_parameters(self):
        Req = self.data['Req']
        Seq = self.data['Seq']
        alpha = self.data['alpha']
        gamma = self.data['gamma']
        sigma = self.data['sigma']
        lambdaa = self.data['lambda']
        delta = self.data['delta']
        mu = self.data['mu']
        
        return (Req, Seq, alpha, gamma, sigma, lambdaa, delta, mu)
    
    # generates a physical model with random parameters, param is the preferred type of parameters drawing
    def generate_random_parameters(self, NR_, NS_, param):
        print('Looking for physical system...')
        model = model_from_random_parameters(NR_, NS_, param)
        repeat = not(model.physical())
        if repeat :
            print("System rejected because it was not physical")
        i=1
        while repeat:
            model = model_from_random_parameters(NR_, NS_, param)
            repeat = not(model.physical())
            i+=1
            if repeat :
                print("System rejected because it was not physical")
        
        print('Physical system found at iteration ', i)
        self.data = model.data
        return 
    
    def set_parameters(self, param):
        self.data = param       
    
    def print_parameters(self):
        print('The model is characterized by the following parameters : ')
        for key in self.data:
            print(key,'=')
            print(self.data[key])
    
    # tests whether the model is physical       
    def physical(self):
        for key in self.data:
            a = self.data[key]
            for i in range(len(np.ravel(a))):
                if np.ravel(a)[i] < 0:
                    return False
        return True
        
    # draws the food network (type ntype) on axis
    def draw_food_network(self, axis, ntype):
        if not('food network' in self.data):
            self.__create_food_network()
            
        colors = dict({'consumption': 'red', 'byproduct release': 'blue'})
        A = self.data['food network']
        
        resources_nodes = [u for (u,d) in A.nodes(data=True) if d['bipartite']==0]
        species_nodes = [u for (u,d) in A.nodes(data=True) if d['bipartite']==1]
    
        resources_sizes = [1000*d['eq_value'] for (u,d) in A.nodes(data=True) if d['bipartite']==0]
        species_sizes = [1000*d['eq_value'] for (u,d) in A.nodes(data=True) if d['bipartite']==1]
    
        relevant_edges = [(u,v) for (u,v,d) in A.edges(data=True) if d['edge_type']==ntype]
        relevant_weights = [3*d['weight'] for (u,v,d) in A.edges(data=True) if d['edge_type']==ntype]
    
        pos = nx.drawing.layout.bipartite_layout(A, species_nodes, align='horizontal')
    
        plt.sca(axis)
        #nx.draw_networkx_labels(A, pos=pos)
        nx.draw_networkx_edges(A, pos=pos, edgelist=relevant_edges, width = relevant_weights, edge_color=colors[ntype])
        nx.draw_networkx_nodes(A, nodelist = resources_nodes, pos=pos, node_size=resources_sizes, node_shape='s')
        nx.draw_networkx_nodes(A, nodelist = species_nodes, pos=pos, node_size=species_sizes, node_shape='h')
        axis.set_title(ntype+' network')
        return
              
    # builds the jacobian from the parameters of the model
    def jacobian(self):
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

        return np.block([[-D, G], [B, Z]])

    def save(self, path):
        np.save(path, self.data)

NR = 3
NS = 3
dissipation = 'random'
mode = ('gamma_construction', 'mainly specialist nested', dissipation)

model = Model()
model.generate_random_parameters(NR, NS, mode)
model.print_parameters()
model.save('test')
model2 = load_model('test.npy')
model2.print_parameters()