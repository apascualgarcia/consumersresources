import matplotlib.pyplot as plt
import numpy as np

folder_path = 'data_output'
filename = "study_feasability_constraints"
savepath = "plots"

fs = 12
error_bar_width = 1.
cap_width = 1.5
markeredgewidth = 2
markersize = 6

data = np.loadtxt(folder_path+'/'+filename+".out", dtype=float)
NR = [int(data[i,1]) for i in range(len(data))]
NS = [int(data[i,2]) for i in range(len(data))]
alpha0 = data[:, 0]
m = [data[i, 3:NR[i]+3] for i in range(len(data))]
d = [data[i, NR[i]+3:] for i in range(len(data))]
print("Data loaded.")

data_meta = [[NR[i], NS[i]] for i in range(len(data))]

# sorted_metaparams[i] is each different set of metaparameters in the file
sorted_metaparams = []
# sorted_m[i] contains all m observed for the set of meta sorted_metaparams
sorted_m = []
sorted_d = []
sorted_alpha0 = []
for i in range(0, len(data_meta)):
    if not(data_meta[i][1] < 0.):
        if data_meta[i] in sorted_metaparams:
            index = sorted_metaparams.index(data_meta[i])
            sorted_m[index].append(m[i].tolist())
            sorted_d[index].append(d[i].tolist())
            sorted_alpha0[index].append(alpha0[i])
        else:
            sorted_metaparams.append(data_meta[i])
            sorted_m.append([m[i].tolist()])
            sorted_d.append([d[i].tolist()])
            sorted_alpha0.append([alpha0[i]])

print("Data sorted.")
for i in range(len(sorted_metaparams)):
    [NR, NS] = sorted_metaparams[i]
    alpha0 = []
    for x in sorted_alpha0[i]:
        if x not in alpha0:
            alpha0.append(x)
    m0 = []
    d0 = []
    error_m0 = []
    error_d0 = []
    for a in alpha0:
        m = np.array([sorted_m[i][j] for j in range(len(sorted_m[i])) if sorted_alpha0[i][j]==a])
        d = np.array([sorted_d[i][j] for j in range(len(sorted_d[i])) if sorted_alpha0[i][j]==a])
        m0.append(np.mean(m))
        d0.append(np.mean(d))
        error_m0.append(np.mean(np.array([np.std(el) for el in m])))
        error_d0.append(np.mean(np.array([np.std(el) for el in d])))

    fig1 = plt.figure('m0 vs a0')
    ax1 = fig1.add_subplot(111)
    ax1.plot(alpha0, m0, '+')
    # ax1.errorbar(alpha0, m0, yerr=error_m0, fmt='+', linestyle='none',
    #              markeredgewidth=markeredgewidth, markersize=markersize,
    #              elinewidth=error_bar_width, capsize=cap_width)
    ax1.set_xlabel(r'$\alpha_0$')
    ax1.set_ylabel(r'$m_0$')
    ax1.set_title(r'Average resource death rate vs. syntrophy (bars indicate average std for one sys)')
    fig1.tight_layout()
    #fig1.savefig(savepath + '/feasability_resource_death_rate.pdf')
    fig1.savefig(savepath + '/feasability_resource_death_rate_no_errorbar.pdf')

    fig2 = plt.figure('d0 vs a0')
    ax2 = fig2.add_subplot(111)
    ax2.plot(alpha0, d0, '+')
    # ax2.errorbar(alpha0, d0, yerr=error_d0, fmt='+', linestyle='none',
    #              markeredgewidth=markeredgewidth, markersize=markersize,
    #              elinewidth=error_bar_width, capsize=cap_width)
    ax2.set_xlabel(r'$\alpha_0$')
    ax2.set_ylabel(r'$d_0$')
    ax2.set_title(r'Average consumer death rate vs. syntrophy (bars indicate average std for one sys)')
    fig2.tight_layout()
    #fig2.savefig(savepath + '/feasability_consumer_death_rate.pdf')
    fig2.savefig(savepath + '/feasability_consumer_death_rate_no_errorbar.pdf')

print('Data plotted and saved')
