import numpy as np
import matplotlib.pyplot as plt
import consumer_resource_data_analysis as cr


folder='data_output'

file='energy_optimal_NIS_matrices'
suptitle=r'NIS'
save_name='plots/energy_NIS'

# file='energy_optimal_LRI_matrices'
# suptitle=r'LRI'
# save_name='plots/energy_LRI'

file='energy_optimal_FC_matrices'
suptitle=r'FC'
save_name='plots/energy_FC'

file='energy_optimal_RS_matrices'
suptitle=r'RS'
save_name='plots/energy_RS'



data=(np.loadtxt(folder+'/'+file+'.out'))
connG=data[:,0]
nestG=data[:,1]
connA=data[:,2]
nestA=data[:,3]
energy=data[:,4]


for i in range(len(connG)):
    connG[i]=cr.closest_element_in_list(connG[i], cr.all_connectance)

for i in range(len(nestG)):
    nestG[i]=cr.closest_element_in_list(nestG[i], cr.all_nestedness)


fig=plt.figure()
ax=fig.add_subplot(111)
for j in range(len(cr.all_nestedness)):
    nest=cr.all_nestedness[j]
    points=[i for i in range(len(connG)) if nestG[i]==nest]
    points=[points[i] for i in np.argsort(connG[points])]
    ax.plot(connG[points], energy[points], color=cr.nest_colours[j], label=r'$\eta_G='+str(nest)+'$')
ax.legend(bbox_to_anchor=(1., 1.))
ax.set_xlabel(r'$\kappa_G$')
ax.set_ylabel(r'$E(A,G)$')
ax.set_title(suptitle)
fig.tight_layout()
fig.savefig(save_name+'_fixed_nest.pdf')

fig=plt.figure()
ax=fig.add_subplot(111)
for j in range(len(cr.all_connectance)):
    conn=cr.all_connectance[j]
    points=[i for i in range(len(nestG)) if connG[i]==conn]
    points=[points[i] for i in np.argsort(nestG[points])]
    ax.plot(nestG[points], energy[points], color=cr.conn_colours[j], label=r'$\kappa_G='+str(conn)+'$')
ax.legend(bbox_to_anchor=(1., 1.))
ax.set_xlabel(r'$\eta_G$')
ax.set_ylabel(r'$E(A,G)$')
ax.set_title(suptitle)
fig.tight_layout()
fig.savefig(save_name+'_fixed_conn.pdf')
