import numpy as np
import copy
import matplotlib.pyplot as plt

gamma_matrix_folder='optimal_matrices/Nr50_Nc25'
gamma_matrix_name='RandTrix_Nr50_Nc25_Nest0.45_Conn0.2768'
# temperature (basically will decide how long the simulation runs)
T = 1
# desired connectance of the alpha matrix
connectance_in = 0.
# maximum number of steps allowed without changing alpha
MaxStepsNotChangedAlpha = 1000
# maximum number of total steps
MaxSteps=7e3
alpha_matrix_name="optimal_alpha_"+gamma_matrix_name


###### DO NOT MODIFY CODE BELOW THIS ########
file_name = gamma_matrix_folder+'/'+gamma_matrix_name+'.txt_corr'
gamma_matrix = np.loadtxt(file_name)

####### FUNCTIONS USED ARE DEFINED HERE ################
def quadratic_form(alpha_):
    to_test = alpha_@gamma_matrix
    return np.linalg.norm(to_test)

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
    qform_=[]
    while not(stop):
        (alpha_,changed)=choose_next_alpha(alpha_,T_)
        connectance_.append(connectance(alpha_))
        qform_.append(quadratic_form(alpha_))
        if not(changed):
            steps_not_changed+=1
        else:
            steps_not_changed=0
        if(steps_not_changed>=MaxStepsNotChangedAlpha):
            stop=True
        steps+=1
        if(steps>=MaxSteps):
            stop=True
    return connectance_, qform_
def random_unit_vector(length_):
    v = np.random.rand(length_)
    return v/np.linalg.norm(v)
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
connect, qdform=MC_algorithm(alpha,T)
np.savetxt(gamma_matrix_folder+'/'+alpha_matrix_name+'.txt', alpha)

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
ax1.set_xlabel('MC Iteration')
ax1.set_ylabel('Connectance')
ax1.plot(connect)
fig1.tight_layout()

fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111)
ax2.set_xlabel('MC Iteration')
ax2.set_ylabel('Quadratic form')
ax2.plot(qdform)
fig2.tight_layout()

fig3 = plt.figure(3)
ax3 = fig3.add_subplot(111)
ax3.set_xlabel('Species')
ax3.set_ylabel('Resources')
ax3.matshow(alpha)
fig3.tight_layout()

plt.show()
