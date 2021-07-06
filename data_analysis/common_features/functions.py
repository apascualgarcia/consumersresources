import numpy as np

# finds the asymptote value of an asymptotic function y, which is assumed to be constant
# or constantly oscillating for its last Npoints
def asymptote(y, Npoints):
    return np.mean(y[-Npoints:])
