#!/usr/local/bin/python3

# import the path for tex to work
import os
os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'
import numpy as np
from scipy.integrate import solve_ivp, odeint
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import matplotlib.style
matplotlib.style.use('default')
np.set_printoptions(precision=2)
np.random.seed(0)


plt.rc('text', usetex=True)
plt.rc('font', family='serif')

NR = 1
NS = 4

R0 = np.ones(NR)*0.4
S0 = np.ones(NS)*3

colors_R=['green', 'olive']
colors_S=['blue', 'red', 'yellow', 'orange']

lambdaa = np.ones(NR)
mu = np.ones(NR)
delta = np.ones(NS)


sigma = np.ones((NS,NR))
gamma = np.ones((NS,NR))
alpha = np.random.uniform(0., 1., (NR,NS))
#alpha = np.ones((NR, NS))
print("alpha = ", alpha)

t0 = 0
tbegin = 90
tend = 100
Npoints = 100

y0 = np.concatenate((R0,S0))
t = np.linspace(tbegin, tend, Npoints)
gammat = np.transpose(gamma)

def norm(a):
    return np.sqrt(np.sum([abs(b)**2 for b in a]))

def abiot_model(t,y):
    R = y[:NR]*1.
    S = y[NR:]*1.
    
    dydt = np.zeros(NR+NS)
    for i in range(NR):
        sum_term = 0.
        for j in range(NS):
            sum_term -= gamma[j][i]*R[i]*S[j]
            sum_term += alpha[i][j]*S[j]
        dydt[i] = lambdaa[i]-mu[i]*R[i]+sum_term
    for i in range(NS):
        sum_term = 0.
        for j in range(NR):
            sum_term += sigma[i][j]*gamma[i][j]*S[i]*R[j]
        dydt[NR+i] = sum_term-delta[i]*S[i]
    return dydt
    
def modified_abiot_model(t,y):
    R = y[:NR]*1.
    S = y[NR:]*1.
    
    dydt = np.zeros(NR+NS)
    
    for i in range(NR):
        sum_term = 0.
        for j in range(NS):
            sum_term -= gamma[j][i]*R[i]*S[j]
            sum_term += alpha[i][j]*S[j]
        dydt[i] = lambdaa[i]-mu[i]*R[i]+sum_term
    
    for i in range(NS):
        sum_term = 0.
        for j in range(NR):
            sum_term += gamma[i][j]*S[i]*R[j] - alpha[j][i]*np.log(R[j])*S[i]
        dydt[NR+i] = sum_term-delta[i]*S[i]
    return dydt
        
    
    return dydt
    
def abiot_model_no_time(y):
    return abiot_model(0, y)
    
def find_equilibrium(yinit, tinit, dt, tol, length=30):
    t0_ = tinit
    tend_ = t0_+length
    sol1 = solve_ivp(abiot_model, (t0_, tend_), yinit, t_eval=[t0_, tend_])
    sol2 = solve_ivp(abiot_model, (tend_, tend_+dt), sol1.y[:,-1], t_eval=[tend_, tend_+dt])
    
    while norm(sol1.y[:-1]-sol2.y[:-1]) < tol:
        t0_  = tend_
        tend_ = t0_+length
        y0_ = sol1.y[:-1]
        sol1 = solve_ivp(abiot_model, (t0_, tend_), y0_, [t0_, tend_])
        sol2 = solve_ivp(abiot_model, (tend, tend+dt), sol1.y[:,-1], [tend_, tend_+dt])
        print('Equilibrium not found after t=', tend_)
    
    print('Equilibrium found after t = ', tend_)
    
    return sol1.y, tend_
    
def alternate_find_equilibrium(estimate):
    sol = fsolve(abiot_model_no_time, estimate)
    return sol
    
def jacobian(y):
    R_ = y[:NR]
    S_ = y[NR:]
    
    M1 = -np.diag(mu+np.dot(gammat, S_))
    M2 = np.dot(-np.diag(R_), gammat)+alpha
    M3 = np.dot(np.diag(S_), np.multiply(sigma, gamma))
    M4 = np.diag(np.dot(np.multiply(sigma, gamma), R_)-delta)
    return np.block([[M1,M2],[M3, M4]])

#yeq, teq = find_equilibrium(y0, t0, dt=20, tol=1e-6)
# sol = solve_ivp(abiot_model, (t0, tend), y0, t_eval = t)
# yeq = alternate_find_equilibrium(sol.y[:,-1])
# jaceq = jacobian(yeq)
# eigenspectrum = np.linalg.eigvals(jaceq)

sol = solve_ivp(modified_abiot_model, (t0, tend), y0, t_eval = t)

for i in range(NR):
    plt.plot(t, sol.y[:][i],'+', color=colors_R[i], markersize=2)
for i in range(NS):
    plt.plot(t, sol.y[:][i+NR],'o', color=colors_S[i], markersize=2)
plt.show()
