import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib as mpl
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erf

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 15}

mpl.rc('font', **font)
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['lines.markersize'] = 8
mpl.rcParams['lines.markeredgewidth'] = 2

save_folder = 'text/figures'
file = 'test'

not_all_conn = False
not_all_nest = False

fs = 12
error_bar_width = 1.
cap_width = 1.5
markeredgewidth = 2
markersize = 6
errorbar = 'no_errorbar'

not_all_conn = True
target_conn = [0.08]

not_all_nest = True
target_nest = [0.1]

def func(x, a,b):
    return erf(-b*(x-a))*0.5+0.5
    #return 1/(1+np.exp(b*(x-a)))

def remove_strings_from_file(filename):
    file = open(filename + '.out', "r")
    metadata = []
    for x in file:
        name = x.replace('./matrices/Nr25_Nc25/RandTrix_Nr', '')
        name = name.replace('_Nc', ' ')
        name = name.replace('_Nest', ' ')
        name = name.replace('_Conn', ' ')
        name = name.replace('.txt', '')
        metadata.append(name)
    file.close()
    f = open(filename + '_filtered.out', 'w')
    for a in metadata:  # python will convert \n to os.linesep
        f.write(a + '\n')
    f.close()

alpha_range = [

]
feasability = [
1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.99,0.99,0.97,0.97,0.99,0.97,0.94,0.9,0.93,0.92,0.89,0.82,0.88,0.86,0.77,0.78,0.71,0.67,0.62,0.57,0.46,0.45,0.51,0.42,0.34,0.31,0.28,0.16,0.23,0.21,0.12,0.18,0.07,0.09,0.07,0.07,0.04,0.04,0.02,0.05,0.02,0.01,0.04,0.01,0.01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
]

popt, pcov = curve_fit(func, alpha_range, feasability)
a,b = popt
print(a)
fitted_curve = [func(x,a,b) for x in alpha_range]
plt.plot(alpha_range, fitted_curve, 'r', label='Fit')
plt.plot(alpha_range, feasability,'b+', label='Data')
plt.legend()
plt.xlabel(r'$\alpha_0$')
plt.ylabel('Probability of feasability')
plt.tight_layout()
plt.savefig(save_folder+'/alpha0_probability_of_feasability.pdf')
plt.show()
