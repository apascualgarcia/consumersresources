import numpy as np
import matplotlib.pyplot as plt

# this script plots the spectrum (and the histogram) of eigenvalues. It saves one file FOR EACH DIFFERENT SET of metaparameters in the file
folder = "/Users/Shared/Cloud/Master/Master_Thesis/Code/consumersresources/TEST_FOLDER/DataOutput"
filename = "eigenvalues_NR25_NS25_0.15_0.1808_eqmag_jac"
savepath = "/Users/Shared/Cloud/Master/Master_Thesis/Code/consumersresources/TEST_FOLDER/Plots"
matrix_path = "Equal_mag_jac"

data = np.loadtxt(folder + "/" + filename + ".out", dtype=complex).tolist()

index_eig_start = 9
# sorted_metaparams[i] is each different set of metaparameters in the file
sorted_metaparams = []
# sorted_eigenvalues[i] contains all eigenvalues observed for the set of meta sorted_metaparams
sorted_eigenvalues = []
for i in range(0, len(data)):
    if not(np.real(data[i][1]) < 0.):
        if data[i][:index_eig_start] in sorted_metaparams:
            index = sorted_metaparams.index(data[i][:index_eig_start])
            sorted_eigenvalues[index].append(data[i][index_eig_start:])
        else:
            sorted_metaparams.append(data[i][:index_eig_start])
            sorted_eigenvalues.append([data[i][index_eig_start:]])
# after detecting all different metaparameters sets, prints spectrum and eigenvalue distribution for all of them
for i in range(len(sorted_metaparams)):

    metaparams = sorted_metaparams[i]
    [gamma0, alpha0, sigma0, p, R0, S0, l0, NR, NS] = np.real(metaparams)

    eigenvalues = np.ravel(sorted_eigenvalues[i])
    color = []
    for j in range(len(eigenvalues)):
        if np.real(eigenvalues[j]) > 1e-16:
            color.append('red')
        else:
            color.append('blue')

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

    fig1 = plt.figure("Eigenspectrum")
    ax1 = fig1.add_subplot(111)
    ax1.set_title(title_string)
    ax1.scatter(np.real(eigenvalues), np.imag(eigenvalues), s=1, c=color)
    ax1.set_xlabel(r'Re($\lambda$)')
    ax1.set_ylabel(r'Im($\lambda$)')

    fig1.savefig(savepath + "/" + matrix_path + "/" + filename +
                 "eigenspectrum_" + save_string + ".png")
    plt.clf()

    fig2 = plt.figure("Real eigenvalue distribution")
    ax2 = fig2.add_subplot(111)
    ax2.set_title(title_string)
    ax2.hist(np.real(eigenvalues), bins=100, density=True)
    fig2.savefig(savepath + "/" + matrix_path + "/" + filename + "real_eigenvalue_distribution_" +
                 save_string + ".pdf")
    fig2.savefig(savepath + "/Real_eigenvalue_spectrum/" + filename + "real_eigenvalue_distribution_" +
                 save_string + ".pdf")
    plt.clf()

    fig3 = plt.figure("Imaginary eigenvalue distribution")
    ax3 = fig3.add_subplot(111)
    ax3.set_title(title_string)
    ax3.hist(np.imag(eigenvalues), bins=100, density=True)
    fig3.savefig(savepath + "/" + matrix_path + "/" + filename + "imag_eigenvalue_distribution_" +
                 save_string + ".pdf")
    fig3.savefig(savepath + "/Imag_eigenvalue_spectrum/" + filename + "imag_eigenvalue_distribution_" +
                 save_string + ".pdf")
    plt.clf()
