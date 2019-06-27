from config import *
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
        dydt[NR+i] = -delta[i]*S[i]+sum_term
    return dydt