#!/usr/local/bin/python3
# import the path for tex to work
import os
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'
import numpy as np
np.set_printoptions(precision=2)
import matplotlib.pyplot as plt

def physical_gamma_sigma(G,S):
    A = np.multiply(S,G)
    return physical_matrix(A)

def physical_matrix(A):
    return np.abs(np.linalg.det(np.dot(np.transpose(A),A))) > 1e-10
    
def create_physical_gamma_sigma(NR, NS):
    gamma = np.random.rand(NS,NR)
    sigma = np.random.rand(NS,NR)
    while not(physical_gamma_sigma(gamma, sigma)):
        gamma = np.random.rand(NS,NR)
        sigma = np.random.rand(NS,NR)
    return gamma, sigma
    
def left_inverse(A):
    return np.dot(np.linalg.inv(np.dot(np.transpose(A), A)), np.transpose(A))
    
    
class Model:    
    def __init__(self, model_param):
        self.NR = model_param['NR']
        self.NS = model_param['NS']
        self.mu = model_param['mu']
        self.alpha = model_param['alpha']
        self.delta = model_param['delta']
        self.gamma = model_param['gamma']
        self.sigma = model_param['sigma']
        self.lambdaa= model_param['lambda']
        self.Req = model_param['Req']
        self.Seq = model_param['Seq']
        
    def physical_model(self):
        Req, Seq = self.equilibria()
        for a in Req:
            if a < 0.:
                return False
        for a in Seq:
            if a < 0:
                return False
        return True
        
        
    def jacobian(self, param):
        R_, S_ = param
        gammat = np.transpose(self.gamma)
        M1 = -np.diag(self.mu+np.dot(gammat, S_))
        M2 = np.dot(-np.diag(R_), gammat)+self.alpha
        M3 = np.dot(np.diag(S_), np.multiply(self.sigma, self.gamma))
        M4 = np.diag(np.dot(np.multiply(self.sigma, self.gamma), R_)-self.delta)
        return np.block([[M1,M2],[M3, M4]])
    
    def compute_equilibria(self):
        if(self.Req=='default' and self.Seq=='default'):
            A = np.multiply(self.sigma, self.gamma)
            AL = left_inverse(A)
            Req = np.dot(AL, self.delta)
        
            B = np.dot(np.diag(Req), np.transpose(self.gamma)) - self.alpha
            BL = left_inverse(B)
            Seq = np.dot(BL, self.lambdaa-np.dot(np.diag(self.mu), Req))
            
            self.Req = Req
            self.Seq = Seq
            
        return
    
    def test_sol(self):
        Req = self.Req
        Seq = self.Seq
        
        res1 = self.lambdaa-np.dot(np.diag(self.mu), Req)-\
        np.dot(np.dot(np.diag(Req), np.transpose(self.gamma)-self.alpha), Seq)
        
        res2 = np.dot(np.dot(np.diag(Seq), np.multiply(self.sigma, self.gamma)), Req)
        res2 = res2-np.dot(np.diag(Seq), self.delta)
        
        return res1, res2

def build_physical_model(NR_, NS_):   
        delta = np.ones(NS)
        mu = np.ones(NR)
        lambdaa = np.ones(NR)
        alpha = np.ones((NR,NS))
        gamma, sigma = create_physical_gamma_sigma()      
        model_params = {
            'NR' : NR_,
            'NS' : NS_,
            'Req' : 'default',
            'Seq' : 'default',
            'mu' : mu,
            'delta' : delta,
            'lambda': lambdaa,
            'alpha' : alpha,
            'gamma' : gamma,
            'sigma' : sigma
        }     
        com = Model(model_params)   
        
        while not(com.physical_model()):
            delta = np.random.uniform(low=0., high=1., size=NS)
            mu = np.random.uniform(low=0., high=1., size=NR)
            lambdaa = np.random.uniform(low=0., high=1., size=NR)
            alpha = np.random.uniform(low=0., high=1., size=(NR,NS))
            gamma, sigma = create_physical_gamma_sigma()
            model_params = {
                'NR' : NR_,
                'NS' : NS_,
                'Req': 'default',
                'Seq' : 'default',
                'mu' : mu,
                'delta' : delta,
                'lambda': lambdaa,
                'alpha' : alpha,
                'gamma' : gamma,
                'sigma' : sigma
            }    
            com = Model(model_params)
        return com
    
def build_model_from_eq(param):
    delta = np.dot(np.multiply(param['sigma'], param['gamma']), param['Req'])
    param['delta'] = delta
    mu = param['lambda']-np.dot(np.dot(np.diag(param['Req']), np.transpose(param['gamma'])), param['Seq'])
    mu = mu + np.dot(param['alpha'], param['Seq'])
    mu = [mu[i]/param['Req'][i] for i in range(len(mu))]
    param['mu'] = mu
    
    com = Model(param)
    return com

def plot_jacobian(Jeq):
    lim = np.max([np.min(Jeq), np.max(Jeq)])    
    plt.matshow(Jeq, cmap='bwr', vmin=-lim, vmax=lim) 
    return

NR = 4
NS = 4

Nsimul = 5000
eigspectrum = []
for i in range(Nsimul):
    Req = np.random.uniform(0.,1.,NR)
    Seq = np.random.uniform(0.,1.,NS)

    model_params = {
        'NR' : NR,
        'NS' : NS,
        'Req': Req,
        'Seq': Seq,
        'lambda' : np.ones(NR),
        'alpha' : np.zeros((NR, NS)),
        'gamma' : np.random.uniform(0.,1., (NS,NR)),
        'sigma' : np.ones((NS,NR)),
    }
    commu = build_model_from_eq(model_params)
    jaceq = commu.jacobian((Req,Seq))
    vals = np.linalg.eigvals(jaceq)
    for j in range(len(vals)):
        eigspectrum.append(vals[j])
print(np.max(eigspectrum))
plt.scatter(np.real(eigspectrum), np.imag(eigspectrum), s=1)
plt.xlim([-2., 0.4])
plt.show()
    