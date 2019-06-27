#!/usr/local/bin/python3
import numpy as np
    
def model_from_random_parameters(NR_, NS_, param):
    (mode, gamma_mode, degree_of_dissipation)=param
    # the idea is to draw beta, Gamma and D such that the stability condition will be easily satisfied
    # and from there reconstruct the rest of the model
    Seq = np.random.uniform(low=0, high=1., size=NS_)
    Req = np.random.uniform(low=0., high=1., size=NR_)
    # take uniform sigma
    sigma = np.random.uniform(low=0.0, high=1., size=(NS_, NR_))
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
            # if repeat :
            #     print("system rejected because Gamma does not fulfill energy dissipation")
    else :
        print('draw_parameters does not know how to implement parameters this way!')

    # find the remaining parameters
    delta = np.dot(np.multiply(sigma, gamma), Req)
    lambdaa = np.dot(np.diag(mu), Req)-np.dot(Gamma, Seq)
    
    model_parameters ={
        'Req': Req,
        'Seq': Seq,
        'alpha': alpha,
        'gamma': gamma,
        'sigma': sigma,
        'lambda': lambdaa,
        'delta': delta,
        'mu': mu
    }
    from microbial_community.microbial_model import Model
    return Model(model_parameters)
    
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
        
