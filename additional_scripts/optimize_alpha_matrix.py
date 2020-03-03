import numpy as np
import copy
import matplotlib.pyplot as plt

gamma_matrix_folder='optimal_matrices/Nr50_Nc25'
gamma_matrix_name='RandTrix_Nr50_Nc25_Nest0.45_Conn0.2768'
# temperature (basically will decide how long the simulation runs)
T = 1
# desired connectance of the alpha matrix
connectance_in = 0.2768
# decides whether or not we allow to release what we do eat
coprophagy = True

###### DO NOT MODIFY CODE BELOW THIS ########
# maximum number of steps allowed without changing alpha
MaxStepsNotChangedAlpha = 1000
# maximum number of total steps
MaxSteps=1e5
# Number of steps needed to record the next set
# of values during the optimization
Nstore=50

alpha_matrix_name="optimal_alpha_"+gamma_matrix_name
file_name = gamma_matrix_folder+'/'+gamma_matrix_name+'.txt'
# first we load the binary alpha matrix
gamma_matrix = np.loadtxt(file_name)
gamma_matrix=np.eye(4)
NS, NR = len(gamma_matrix), len(gamma_matrix[0])
# then we define the unitary vector with which the quadratic
# form will be computed
unit = np.ones(NR)


####### FUNCTIONS USED ARE DEFINED HERE ################
def quadratic_form(alpha_):
    to_test = alpha_@gamma_matrix
    return np.transpose(unit)@to_test@unit

def probability_density(alpha_, T_):
    return np.exp(-quadratic_form(alpha_)/T_)

def proposed_new_alpha(alpha_, steps_):
    return proposed_new_alpha_flip(alpha_, steps_)

def choose_next_alpha(alpha_, T_, steps_, fails_):
    new_alpha_=proposed_new_alpha(alpha_, steps_)
    proba_ratio_=probability_density(new_alpha_, T_)/probability_density(alpha_, T_)
    if proba_ratio_ > 1:
        # in that case the move is optimal and is hence accepted
        alpha_ = copy.deepcopy(new_alpha_)
        fails_ = 0
    else:
        fails_+=1
        # we may also accept a suboptimal move
        if(np.random.uniform() < proba_ratio_):
            alpha_=copy.deepcopy(new_alpha_)

    # otherwise we do not accept the move and keep the previous alpha
    return

def MC_algorithm(alpha_):
    stop_=False
    steps_=0
    fails_=0
    while not(stop_):
        # CAREFUL DID NOT IMPLEMENT THE FAIL ARGUMENT YET
        choose_next_alpha(alpha_,T, steps_, fails_)
        if(steps_%100==0):
            print(quadratic_form(alpha_))
            print(alpha_)
        steps_+=1
        if(steps_>=MaxSteps):
            stop_ = True
    return

def random_unit_vector(length_):
    v = np.random.rand(length_)
    return v/np.linalg.norm(v)

def create_alpha(connectance_in_, gamma_):
    # we want the connectance of alpha to be connectance_
    alpha_ = np.zeros((NR, NS))
    for mu in range(NR):
        for i in range(NS):
            if(np.random.uniform()<connectance_in_):
                alpha_[mu][i]=1
    return alpha_

def proposed_new_alpha_Alberto(alpha_, steps_):
    # we alternate between modifying rows and columns
    if steps_%2==0:
        new_alpha_ = modify_row(alpha_, gamma_matrix, coprophagy)
    else:
        new_alpha_ = modify_column(alpha_, gamma_matrix, coprophagy)
    return new_alpha_

def modify_row(alpha_, gamma_, coprophagy_):
    # we first choose a random row (i.e resource) that is not empty and not
    # completely filled
    new_alpha_ = copy.deepcopy(alpha_)
    sum_row=np.sum(new_alpha_, axis=1)
    mu=np.random.randint(low=0, high=NR)
    while(sum_row[mu]==0 or sum_row[mu]==NR):
        mu=np.random.randint(low=0, high=NR)

    # then we find the indices of the columns (i.e. consumers) in that row
    # (i.e. resource) which are zeros and which are ones
    zero_els=[]
    one_els=[]
    for i in range(NS):
        if new_alpha_[mu][i]==0:
            zero_els.append(i)
        elif new_alpha_[mu][i]==1:
            one_els.append(i)

    # we choose one index which has one and another which has zero
    zero_index=np.random.choice(zero_els)
    one_index=np.random.choice(one_els)

    # we check if there is another non empty value in the column (i.e.
    # check if that species releases to something else)
    conditionRel=False
    for nu in range(NR):
        if new_alpha_[nu][one_index]==1:
            conditionRel=True

    # we check if with the chosen zero index, there is coprophagy
    conditionCopr=False
    if(gamma_[zero_index][mu]!=1):
        conditionCopr=True
    elif coprophagy_==True:
        conditionCopr=True

    # if both conditions are fulfilled we can swap the two elements
    if (conditionCopr and conditionRel):
        new_alpha_[mu][zero_index]=1
        new_alpha_[mu][one_index]=0

    return new_alpha_

def modify_column(alpha_, gamma_, coprophagy_):
    # we first choose a random column (i.e consumer) that is not empty and not completely filled
    new_alpha_ = copy.deepcopy(alpha_)
    sum_column = np.sum(new_alpha_, axis=0)
    k = np.random.randint(low=0, high=NS)
    while((sum_column[k]==0) or (sum_column[k])==NS):
        k = np.random.randint(low=0, high=NS)

    # then we find the indices of the rows in that column which are zeros or
    # which are one
    zero_els = []
    one_els = []
    for mu in range(NR):
        if new_alpha_[mu][k]==0:
            zero_els.append(mu)
        elif new_alpha_[mu][k]==1:
            one_els.append(mu)

    # we then choose randomly one of the zero elements and one of the one elements
    zero_index = np.random.choice(zero_els)
    one_index = np.random.choice(one_els)

    # we check if there is another non empty value in the row with the one, i.e.
    # i.e. if the resource is being released by another species
    conditionRel=False
    for i in range(NS):
        if(alpha[one_index][i]==1):
            conditionRel=True

    # we check if with the chosen zero index, there is coprophagy
    conditionCopr=False
    if(gamma_[k][zero_index]!=1):
        conditionCopr=True
    elif coprophagy_==True:
        conditionCopr=True

    # if both conditions are fulfilled we can swap the two elements
    if (conditionCopr and conditionRel):
        new_alpha_[zero_index][k]=1
        new_alpha_[one_index][k]=0

    return new_alpha_

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

def proposed_new_alpha_flip(alpha_, steps_):
    new_alpha_=copy.deepcopy(alpha_)
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
# MAIN IS HERE
# first we create an alpha matrix
alpha=create_alpha(connectance_in, gamma_matrix)

# then we apply the Monte Carlo algorithm to it
MC_algorithm(alpha)
