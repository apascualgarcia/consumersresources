import numpy as np
import matplotlib.pyplot as plt



def f(l, alpha, beta):
    sum=0.
    for i in range(len(alpha)):
        sum+=np.log(l**2+alpha[i]*l-beta[i])
    return sum-1
    

