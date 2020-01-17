import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm

folder_path = "/Users/Shared/Cloud/Master/Master_Thesis/Code/consumersresources/TEST_FOLDER/DataOutput"
jacob_filename = "jacobian_NR25_NS50_0.4_0.3848_equal_jac"
savepath = "/Users/Shared/Cloud/Master/Master_Thesis/Code/consumersresources/TEST_FOLDER/Plots"
matrix_path = "Equal_mag_jac"

data_meta = np.loadtxt(folder_path + "/" +
                       jacob_filename + ".out.meta").tolist()
data_jac = np.loadtxt(folder_path + "/" + jacob_filename + ".out").tolist()

print("Data loaded.")

# sorted_metaparams[i] is each different set of metaparameters in the file
sorted_metaparams = []
# sorted_eigenvalues[i] contains all eigenvalues observed for the set of meta sorted_metaparams
sorted_jacobian = []
for i in range(0, len(data_meta)):
    if not(data_meta[i][1] < 0.):
        NR = int(data_meta[i][7])
        NS = int(data_meta[i][8])
        jac = data_jac[i * (NR + NS):(i + 1) * (NR + NS)][:]
        if data_meta[i] in sorted_metaparams:
            index = sorted_metaparams.index(data_meta[i])
            sorted_jacobian[index].append(jac)
        else:
            sorted_metaparams.append(data_meta[i])
            sorted_jacobian.append([jac])
print("Data sorted.")

minval = np.max(np.array(sorted_jacobian))
maxval = np.min(np.array(sorted_jacobian))
mean_matrices = []

# after detecting all different metaparameters sets, compute mean jacobian matrix for all of them
for i in range(len(sorted_metaparams)):
    jacobians = sorted_jacobian[i]
    Nsimul = len(jacobians)
    length = len(jacobians[0])

    if(length != len(jacobians[0][0])):
        print("Careful, matrices ill formed")

    mean_matrix = np.zeros(shape=(length, length))
    for j in range(length):
        for k in range(length):
            for l in range(Nsimul):
                mean_matrix[j][k] += jacobians[l][j][k] / Nsimul
    if(np.min(mean_matrix) < minval):
        minval = np.min(mean_matrix)
    if(np.max(mean_matrix) > maxval):
        maxval = np.max(mean_matrix)
    mean_matrix = np.ma.masked_where(np.abs(mean_matrix) < 1e-15, mean_matrix)
    mean_matrices.append(mean_matrix)
print("Average matrices computed")

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

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(title_string)
    pos = ax.imshow(mean_matrices[i], cmap=plt.cm.get_cmap(
        'jet'), norm=SymLogNorm(linthresh=0.01, vmin=minval, vmax=maxval))
    fig.colorbar(pos)
    fig.savefig(savepath + "/Jacobians/" +
                jacob_filename + "_" + save_string + ".png")
    fig.savefig(savepath + "/" + matrix_path + "/" +
                jacob_filename + "_" + save_string + ".png")
    plt.clf()
print("Average matrices plotted and saved.")
