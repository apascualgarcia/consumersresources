import matplotlib.pyplot as plt
import numpy as np

save_folder='data_output'
filename='critical_delta_convergence_study_conv_configuration_comparison_NR25_NS25_s05_a0_random_matrix_NR25_NS25_S0=0.1_gamma0=0.5'
matrices_folder='./matrices/Nr25_Nc25'


def remove_strings_from_file(filename):
    file = open(filename + '.out', "r")
    metadata = []
    for x in file:
        name = x.replace(matrices_folder + '/RandTrix_Nr', '')
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

remove_strings_from_file(save_folder+'/'+filename)
data = np.loadtxt(save_folder+'/'+filename+'_filtered.out')

NR=data[0,0]
NS=data[0,1]
nestedness=data[0,2]
connectance=data[0,3]
convergence_threshold=data[:,4]
critical_delta=data[:,5]
critical_delta_error=data[:,6]

fig1=plt.figure(1)
ax1=fig1.add_subplot(111)
title=r'$N_R='+str(NR)+', N_S='+str(NS)+', \eta='+str(nestedness)+', \kappa='+str(connectance)+'$'
ax1.plot(convergence_threshold, critical_delta)
ax1.set_title(title)
ax1.set_xscale('log')
ax1.set_xlabel(r'$\epsilon_{thr}$')
ax1.set_ylabel(r'$\Delta^*$')
fig1.tight_layout()

fig1.savefig('./plots/'+filename+'.pdf')

plt.show()
