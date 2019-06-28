import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from microbial_community.parameters_draw import model_from_random_parameters
from scipy.integrate import solve_ivp
from matplotlib.lines import Line2D

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
        # if repeat :
        #     print("System rejected because it was not physical")
        i=1
        while repeat:
            model = model_from_random_parameters(NR_, NS_, param)
            repeat = not(model.physical())
            i+=1
            # if repeat :
            #         print("System rejected because it was not physical")
        
        print('Physical system found at iteration ', i)
        self.data = model.data
        self.data['NR'] = NR_
        self.data['NS'] = NS_
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
        
    #draws the food network (type ntype) on axis
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

        #nx.draw_networkx_labels(A, pos=pos)
        nx.draw_networkx_edges(A, pos=pos, edgelist=relevant_edges, width = relevant_weights, edge_color=colors[ntype], ax = axis)
        nx.draw_networkx_nodes(A, nodelist = resources_nodes, pos=pos, node_size=resources_sizes, node_shape='s', ax = axis)
        nx.draw_networkx_nodes(A, nodelist = species_nodes, pos=pos, node_size=species_sizes, node_shape='h', ax = axis)
        axis.set_title(ntype+' network')
        axis.axis('on')
        axis.set_xticklabels([])
        axis.set_yticklabels([])
        axis.set_xticks([])
        axis.set_yticks([])
        return

    # builds the jacobian from the parameters of the model
    def jacobian(self):
        (Req, Seq, alpha, gamma, sigma, lambdaa, delta, mu) = self.retrieve_model_parameters()
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
            
        J = np.block([[-D, G], [B, Z]])
        
        self.data['Jeq'] = J
        return J

    # returns the solution to the time evolution of the system 
    # taken into account we started at time t=0 with the densities Rinit, Sinit
    def time_evolution(self, Rinit, Sinit, t):
        y0 = np.zeros(len(Rinit)+len(Sinit))
        
        for i in range(len(Rinit)):
            y0[i] = Rinit[i]
            
        for i in range(len(Sinit)):
            y0[i+len(Rinit)]=Sinit[i]
        
        y0 = np.array(y0)
        sol = solve_ivp(self.__system_evolve, t_span=(0, t[-1]), y0=y0, t_eval=t)
        
        return (sol.y[:self.data['NR']], sol.y[self.data['NR']:], sol.t) 
        
    def __system_evolve(self,t,env):
        (NR, lambdaa, mu, gamma, alpha, sigma, delta) = (self.data['NR'], self.data['lambda'], self.data['mu'], self.data['gamma'],\
        self.data['alpha'], self.data['sigma'], self.data['delta'])
        
        R = env[:NR]
        S = env[NR:]
    
        denvdt = np.zeros(len(R)+len(S))
            
        for i in range(len(R)):
            denvdt[i] = lambdaa[i]-mu[i]*R[i]-R[i]*sum([gamma[j,i]*S[j] for j in range(len(S))])+sum([alpha[i,j]*S[j] for j in range(len(S))])
            
        for i in range(len(S)):
            denvdt[i+len(R)] = (sum([sigma[i,j]*gamma[i,j]*R[j] for j in range(len(R))])-delta[i])*S[i]
        
        return denvdt
        
    def save(self, path):
        np.save(path, self.data)
        return
        
    def load(self, filepath):
        data = np.load(filepath+'.npy', allow_pickle=True)
        self.data = data[()]
        return
        
    # plots the spectrum of the jacobian
    def plot_eigenspectrum(self, axis):
        
        eigenvalues = np.linalg.eigvals(self.jacobian())
        axis.scatter(np.real(eigenvalues), np.imag(eigenvalues), color = 'red')
        axis.plot([0., 0.], [np.max(np.imag(eigenvalues)), np.min(np.imag(eigenvalues))], color = 'grey')
        
        lambda_max = np.max(np.real(eigenvalues))
        titre = 'Spectrum of the jacobian at equilibrium '
        if lambda_max > 0:
            titre += '(unstable)'
        else:
            titre += '(stable)'
            
        axis.set_title(titre)
        axis.set_xlabel(r'Re($\lambda$)')
        axis.set_ylabel(r'Im($\lambda$)')
        
        return
        
        
        

# takes the output of model.time_involution and plots the species and resources joint time evolution on the given axis
def plot_time_evolution(axis, sol):
    (R, S, t) = sol
        
    for i in range(len(R)):
        axis.plot(t, R[i], color = 'green')
    for i in range(len(S)):
        axis.plot(t, S[i], color = 'blue')
        
    custom_lines = [Line2D([0], [0], color='green', lw=4),\
    Line2D([0], [0], color='blue', lw=4)]
    axis.axis('on')
    axis.legend(custom_lines, ['Resources', 'Species'])
    axis.set_title('Time evolution of the ecosystem')
    
    return 
    
    