import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch

# folder='data_output'
# file='optimal_LRI_for_NR25_NS25'
# Nest=['0.599869',  '0.104794']
# Conn=['0.316800',  '0.129600']
Gamma_matrix=['optimal_matrices/consumption/Nr25_Nc25/RandTrix_Nr25_Nc25_Nest0.6_Conn0.3168.txt',
              'optimal_matrices/consumption/Nr25_Nc25/RandTrix_Nr25_Nc25_Nest0.1_Conn0.1296.txt']
Alpha_matrix=['optimal_matrices/syntrophy/optimal_LRI_Nr25_Nc25/RandTrix_Nr25_Nc25_Nest0.6_Conn0.3168_optimal_alpha.txt',
              'optimal_matrices/syntrophy/optimal_LRI_Nr25_Nc25/RandTrix_Nr25_Nc25_Nest0.1_Conn0.1296_optimal_alpha.txt']

mat_g = []
mat_a = []
for i in range(2):
    mat_g.append(np.loadtxt(Gamma_matrix[i]))
    mat_a.append(np.loadtxt(Alpha_matrix[i]))

fig, axs= plt.subplots(2,2)
for i in range(2):
    axs[i,0].imshow(mat_g[i], extent=[1,len(mat_g[i]), len(mat_g[i][0]),1], cmap='Greys')
    axs[i,0].set_title(r'$G_'+str(i+1)+'$')
    axs[i,0].set_xlabel(r'$N_R$')
    axs[i,0].set_ylabel(r'$N_S$')

    axs[i,1].imshow(mat_a[i], extent=[1, len(mat_a[i]), len(mat_a[i][0]), 1], cmap='Greys')
    axs[i,1].set_title(r'$A_'+str(i+1)+'$')
    axs[i,1].set_ylabel(r'$N_R$')
    axs[i,1].set_xlabel(r'$N_S$')
    #axs[i,1].set_title(r'$\kappa_G='+str(round(float(Conn[i]),2))+', \eta_G='+str(round(float(Nest[i]),2))+'$')
    xy1=(30,12.5)
    xy2 = (-15, 12.5)
    con = ConnectionPatch(xyA=xy1, xyB=xy2, coordsA=axs[i,0].transData, coordsB=axs[i,1].transData,
                       color="black", lw=2, arrowstyle="->")
    axs[i,0].add_artist(con)
fig.subplots_adjust(hspace=0.01)
fig.tight_layout()
fig.savefig('plots/typical_optimal_LRI_matrix.pdf')
plt.close()


fig2, ax2 = plt.subplots(2,2)
vmin=0
vmax=18
for i in range(2):
    AG = mat_a[i]@mat_g[i]
    GTGAG= np.transpose(mat_g[i])@mat_g[i]-0.16*mat_a[i]@mat_g[i]
    ax2[i,0].imshow(AG, cmap='jet',vmin=0, vmax=18, extent=[1, len(AG), len(AG[0]), 1])
    im=ax2[i,1].imshow(GTGAG, cmap='jet', vmin=0, vmax=18, extent=[1, len(AG), len(AG[0]), 1])

    ax2[i,0].set_title(r'$AG$')
    ax2[i,1].set_title(r'$G^TG-\alpha_0/(\gamma_0 R_0)AG$')

    ax2[i,0].set_xlabel(r'$N_R$')
    ax2[i,0].set_ylabel(r'$N_R$')
    ax2[i,1].set_xlabel(r'$N_R$')
    ax2[i,1].set_ylabel(r'$N_R$')
fig2.subplots_adjust(bottom=0.7, hspace=0.1)
cbar_ax = fig2.add_axes([0.15, 0.05, 0.7, 0.02])
fig2.colorbar(im, cax=cbar_ax, orientation='horizontal')
fig2.tight_layout()
plt.show()
