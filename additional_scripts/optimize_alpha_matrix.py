import numpy as np
import copy
import matplotlib.pyplot as plt

gamma_matrix_folder='optimal_matrices/Nr25_Nc25'
gamma_matrix_name='RandTrix_Nr25_Nc25_Nest0.2_Conn0.1312'
# temperature (basically will decide how long the simulation runs)
T = 1.
# desired connectance of the alpha matrix
connectance_in = 1.
# maximum number of steps allowed without changing alpha
MaxStepsNotChangedAlpha = 1000
# maximum number of total steps
MaxSteps=1e6

file_name = gamma_matrix_folder+'/'+gamma_matrix_name+'.txt'
#gamma_matrix = np.loadtxt(file_name)
gamma_matrix = np.array([
[1,1,1,1],
[1,1,1,0],
[1,1,0,0],
[1,0,0,0]
])

####### FUNCTIONS USED ARE DEFINED HERE ################
def quadratic_form(alpha_):
    to_test = alpha_@gamma_matrix
    to_test = 0.5*(to_test+np.transpose(to_test))
    return np.max(np.linalg.eigvalsh(to_test))
def probability_density(alpha_, T_):
    return np.exp(-quadratic_form(alpha_)/T_)
def proposed_new_alpha(alpha_):
    return proposed_new_alpha_flip(alpha_)
def choose_next_alpha(alpha_, T_):
    new_alpha_=proposed_new_alpha(alpha_)
    proba_ratio_=probability_density(new_alpha_, T_)/probability_density(alpha_, T_)
    if proba_ratio_ > 1:
        return (new_alpha_, True)
    else:
        if(np.random.uniform() < proba_ratio_):
            return (new_alpha_, True)
    return (alpha_, False)
def MC_algorithm(alpha_, T_):
    stop=False
    steps=0
    steps_not_changed=0
    connectance_=[]
    while not(stop):
        (alpha_,changed)=choose_next_alpha(alpha_,T_)
        if not(changed):
            steps_not_changed+=1
        else:
            steps_not_changed=0
        if(steps_not_changed>=MaxStepsNotChangedAlpha):
            stop=True
        steps+=1
        if(steps>=MaxSteps):
            stop=True
        connectance_.append(connectance(alpha_))
    return
def create_alpha(connectance_in, gamma_):
    # we want the connectance of alpha to be connectance_
    rows = len(gamma_[0])
    cols = len(gamma_)
    alpha_ = np.zeros((rows, cols))
    for mu in range(rows):
        for i in range(cols):
            if(np.random.uniform()<connectance_in):
                alpha_[mu][i]=1
    return alpha_
def proposed_new_alpha_exchange(alpha_):
    # either change two rows or two columns
    new_alpha_ = alpha_
    to_swap=[0,0]

    if(np.random.uniform()<0.5):
        rows = len(alpha_)
        to_swap = np.random.randint(0, rows, size=2)
        while to_swap[0]==to_swap[1]:
            to_swap = np.random.randint(0, rows, size=2)
        exchange_rows(to_swap[0], to_swap[1], new_alpha_)
    else:
        cols = len(alpha_[0])
        to_swap = np.random.randint(0, cols, size=2)
        while to_swap[0]==to_swap[1]:
            to_swap = np.random.randint(0, cols, size=2)
        exchange_cols(to_swap[0], to_swap[1], new_alpha_)
    return new_alpha_
def proposed_new_alpha_flip(alpha_):
    new_alpha_=alpha_
    # simply take a random element and flip it
    rows = len(alpha_)
    cols = len(alpha_[0])
    chosen_row, chosen_col = np.random.randint(0, rows), np.random.randint(0, cols)
    new_alpha_[chosen_row][chosen_col]=1-new_alpha_[chosen_row][chosen_col]
    return new_alpha_
def connectance(alpha_):
    dimension_ = len(alpha_)*len(alpha_[0])
    connectance_=0
    for line in alpha:
        for el in line:
            if el>0:
                connectance_+=1
    return connectance_/dimension_
def exchange_rows(i_, j_, mat_):
    mat_[[i_,j_]]=mat_[[j_,i_]]
    return
def exchange_cols(i_, j_, mat_):
    exchange_rows(i_, j_, np.transpose(mat_))
    return

####################################################


############ MAIN IS HERE ##########################
alpha = create_alpha(connectance_in, gamma_matrix)
connectance=MC_algorithm(alpha,T)
print(connectance)
