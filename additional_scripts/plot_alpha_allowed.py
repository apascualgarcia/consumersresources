import matplotlib.pyplot as plt
import numpy as np


def alpha_threshold(eps, s0, g0, R0, S0, l0, NR, NS):
    a = (1 - eps)**2 / (1 + eps) * (1 - (1 + eps) * s0) * g0 * R0 / NR
    b = (1 - eps)**3 / (1 + eps) * s0 * g0 * R0 / NR
    c = (NS * g0 * R0 * S0 * (1 + eps)**3 -
         l0 * (1 - eps)) / (S0 * (1 - eps)**2)
    return [c, np.min(np.array([a, b]))]


def alpha_allowed(alpha, eps, s0, g0, R0, S0, l0, NR, NS):
    min_val, max_val = alpha_threshold(eps, s0, g0, R0, S0, l0, NR, NS)
    if alpha > min_val and alpha < max_val:
        return 1.
    return 0.


def C_range(eps, s0, g0, R0, S0, l0, NR, NS):
    a = R0 / (s0 * S0)
    b = 1 / s0 * (R0 / S0 * (1 - NS * (1 + eps)**3 /
                             ((1 - eps)**2)) + l0 / (g0 * S0**2) * 1 / 1 - eps)
    c = 0
    if s0 <= 0.5:
        c = R0 / S0 * (1 / s0 - (1 - eps)**3 / ((1 + eps) * NR))
    else:
        c = R0 / S0 * (1 / s0 * (1 - (1 - eps)**2 /
                                 (NR * (1 + eps))) + (1 - eps)**2 / NR)
    return [c, np.min([a, b])]


def C_value(eps, s0, g0, R0, S0, l0, NR, NS):
    a, b = C_range(eps, s0, g0, R0, S0, l0, NR, NS)
    return (a + b) * 0.5


eps = 0.1
g0 = 1.
s0 = 0.5
R0 = 300
S0 = 1.0
l0 = 11092

NR = 25
NS = 25

print(alpha_threshold(eps, s0, g0, R0, S0, l0, NR, NS))
