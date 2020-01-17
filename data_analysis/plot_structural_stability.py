import matplotlib.pyplot as plt
import numpy as np

filename = "DataOutput/structural_stability_Nr25_Nc25_Nest0.15_Conn0.1808.out"
data = np.loadtxt(filename).tolist()

# sorted_metaparams[i] is each different set of metaparameters in the file
metaparams = []
# sorted_m[i] contains all m observed for the set of meta sorted_metaparams
time_eq = []
extinctions = []

data_meta = [x[0:3] for x in data]

for i in range(0, len(data_meta)):
    if data_meta[i] in metaparams:
        index = metaparams.index(data_meta[i])
        time_eq[index].append(data[i][3])
        extinctions[index].append(data[i][4])
    else:
        metaparams.append(data_meta[i])
        time_eq.append([data[i][3]])
        extinctions.append([data[i][4]])
print("Data sorted.")

average_extinctions = [
    np.mean(extinctions[i]) for i in range(len(extinctions))]
average_time_eq = [np.mean(time_eq[i])
                   for i in range(len(time_eq))]

sorted_metaparams = []
sorted_delta = []
sorted_time_eq = []
sorted_extinctions = []

for i in range(0, len(metaparams)):
    mp = metaparams[i][0:2]
    if mp in sorted_metaparams:
        index = sorted_metaparams.index(mp)
        sorted_delta[index].append(metaparams[i][2])
        sorted_time_eq[index].append(average_time_eq[i])
        sorted_extinctions[index].append(average_extinctions[i])
    else:
        sorted_metaparams.append(mp)
        sorted_delta.append([metaparams[i][2]])
        sorted_time_eq.append([average_time_eq[i]])
        sorted_extinctions.append([average_extinctions[i]])

for i in range(len(sorted_metaparams)):
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(sorted_delta[i], sorted_extinctions[i], '+-')
    ax1.set_xlabel(r'$\Delta$')
    ax1.set_ylabel("Number of extinctions")

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(sorted_delta[i], sorted_time_eq[i], '+-')
    ax2.set_xlabel(r'$\Delta$')
    ax2.set_ylabel("Time to reach equilibrium")

    plt.show()
    plt.clf()
