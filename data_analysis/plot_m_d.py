import matplotlib.pyplot as plt
import numpy as np

folder_path = "/Users/Shared/Cloud/Master/Master_Thesis/Code/consumersresources/TEST_FOLDER/DataOutput"
filename = "death_rates_NR25_NS50_0.4_0.3848_equal_jac"
savepath = "/Users/Shared/Cloud/Master/Master_Thesis/Code/consumersresources/TEST_FOLDER/Plots"
matrix_path = "Equal_mag_jac"

m = np.loadtxt(folder_path + "/" + filename + ".out.resources").tolist()
d = np.loadtxt(folder_path + "/" + filename + ".out.consumers").tolist()

data_meta = np.loadtxt(folder_path + "/" + filename + ".out.meta").tolist()
print("Data loaded.")


# sorted_metaparams[i] is each different set of metaparameters in the file
sorted_metaparams = []
# sorted_m[i] contains all m observed for the set of meta sorted_metaparams
sorted_m = []
sorted_d = []
for i in range(0, len(data_meta)):
    if not(data_meta[i][1] < 0.):
        if data_meta[i] in sorted_metaparams:
            index = sorted_metaparams.index(data_meta[i])
            sorted_m[index].append(m[i])
            sorted_d[index].append(d[i])
        else:
            sorted_metaparams.append(data_meta[i])
            sorted_m.append([m[i]])
            sorted_d.append([d[i]])
print("Data sorted.")

for i in range(len(sorted_metaparams)):
    metaparams = sorted_metaparams[i]
    [gamma0, alpha0, sigma0, p, R0, S0, l0, NR, NS] = np.real(metaparams)
    title_string = r"$\gamma_0=$" + str(gamma0)
    title_string += r", $\sigma_0=$" + str(sigma0)
    title_string += r", $\alpha_0=$" + str(alpha0)
    title_string += r", $N_R=$" + str(NR)
    title_string += r", $N_S=$" + str(NS)

    save_string = title_string.replace('_', '')
    save_string = save_string.replace(" ", "_")
    save_string = save_string.replace('$', '')
    save_string = save_string.replace(',', '')
    save_string = save_string.replace('\\', '')

    s_m = np.ravel(np.array(sorted_m[i]))
    s_d = np.ravel(np.array(sorted_d[i]))

    fig1 = plt.figure("Resources death rate distribution")
    ax1 = fig1.add_subplot(111)
    ax1.hist(s_m, bins=100, density=True)
    ax1.set_xlabel(r'$m_\mu$')
    ax1.set_ylabel("Probability")
    ax1.set_title(title_string)
    fig1.savefig(savepath + "/Death_rates_resources/" + filename +
                 "_resources_death_rate" + save_string + ".pdf")
    fig1.savefig(savepath + "/" + matrix_path + "/" + filename +
                 "_resources_death_rate" + save_string + ".pdf")
    plt.clf()

    fig2 = plt.figure("Consumers death rate distribution")
    ax2 = fig2.add_subplot(111)
    ax2.hist(s_d, bins=100, density=True)
    ax2.set_title(title_string)
    ax2.set_xlabel(r'$d_i$')
    ax2.set_ylabel("Probability")
    fig2.savefig(savepath + "/" + matrix_path + "/" + filename +
                 "_consumers_death_rate" + save_string + ".pdf")
    plt.clf()
print("Death rate distribution plotted and saved.")
