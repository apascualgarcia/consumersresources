import matplotlib.pyplot as plt
import numpy as np
import consumer_resource_data_analysis
from scipy.optimize import curve_fit

filename = "data_output/typical_probability_structural_stability_curve.out"
data = np.loadtxt(filename)

delta = data[0::2]
proba = data[1::2]

def sigmoidal_fit(x, x0, s):
    return 1./(1.+np.exp(-s*(x-x0)))



fig = plt.figure()
ax = fig.add_subplot(111)
popt, pcov = curve_fit(sigmoidal_fit, delta, proba, maxfev=9000000)
fitted_curve = [sigmoidal_fit(x, *popt) for x in delta]

ax.set_xlabel(r'$\Delta_S$')
ax.set_ylabel(r'$P_E(\Delta_S)$')
ax.set_xlim(min(delta)-0.002, max(delta)+0.002)
ax.set_ylim(-0.05, 1.05)

xdata, ydata=popt[0],-0.05

bbox = dict(boxstyle="square", fc="0.99")
arrowprops = dict(
    arrowstyle = "->",
    connectionstyle = "arc3, rad=0",
    relpos=(0,0))

offset = 40
ax.annotate('$\Delta_S^* = %.3f $'%xdata,
            (xdata, ydata), xytext=(offset, offset), textcoords='offset points',
            bbox=bbox, arrowprops=arrowprops)
ax.plot(delta, fitted_curve, linestyle='solid', marker='', color='red', label='Sigmoidal fit')
ax.plot(delta, proba, linestyle='None', color='blue', label='Numerical data')

ax.plot(np.linspace(0, popt[0],10), [0.5 for x in range(10)], linestyle='dashed', marker='', color='black')
ax.plot([popt[0] for x in range(10)], np.linspace(-0.05, 0.5, 10), linestyle='dashed', marker='', color='black')

ax.set_yticks(np.linspace(0,1,5))
labels =  np.round(np.linspace(0,1,5), decimals=2)
ax.set_yticklabels([r'$'+str(a)+'$' for a in labels])

ax.legend()

fig.tight_layout()
fig.savefig("plots/typical_probability_structural_stability_curve.pdf")

plt.show()
