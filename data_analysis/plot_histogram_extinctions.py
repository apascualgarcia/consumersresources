import numpy as np
import matplotlib.pyplot as plt

# nest = 0.1
# conn = 0.0832
#
# nest = 0.25
# conn = 0.2336
#
# nest = 0.55
# conn = 0.2784
#
nest = 0.55
conn = 0.4176

filename = 'histogram_extinctions_nest' + str(nest) + '_conn' + str(conn)
data = np.loadtxt('data_output/' + filename + ".out")

extinct = data[:, 0]
time_eq = data[:, 1]


mean_extinct = round(np.mean(extinct), 2)
mean_time_eq = round(np.mean(time_eq), 2)

title_string = 'Nestedness ' + \
    str(nest) + ' Connectance ' + str(conn) + \
    ' (' + str(len(extinct)) + ' runs)'
nbins = np.linspace(min(extinct) - 0.5, max(extinct) +
                    0.5, max(extinct) - min(extinct) + 2)
fig1 = plt.figure('Histogram extinctions')
ax1 = fig1.add_subplot(111)
ax1.hist(extinct, nbins, align='mid', density=1)
ax1.set_title(title_string + '(mean is ' + str(mean_extinct) + ')')
ax1.set_xlabel(r'Number of extinctions at $\Delta = \Delta^*$')
ax1.set_ylabel(r'Probability')
ax1.set_xticks(range(int(min(extinct)), int(max(extinct) + 1)))
fig1.savefig('plots/' + filename + '.pdf')


fig2 = plt.figure('Histogram time to reach equilibrium')
ax2 = fig2.add_subplot(111)
bins = np.logspace(np.log10(min(time_eq)), np.log10(max(time_eq)), 100)
# bins = 100
ax2.hist(time_eq, bins, align='mid', density=1)
ax2.set_title(title_string + '(mean is ' + str(mean_time_eq) + ')')
ax2.set_xlabel(
    r'Time to reach equilibrium after perturbation at $\Delta = \Delta^*$')
ax2.set_ylabel(r'Probability')
ax2.set_xscale('log')
ax2.set_yscale('log', nonposy='clip')
fig2.savefig('plots/' + 'histogram_teq_nest' +
             str(nest) + '_conn' + str(conn) + '.pdf')
plt.clf()
