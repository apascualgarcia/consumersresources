import numpy as np
from consumer_resource_data_analysis import all_nestedness, all_connectance, alpha_mode_colours
import consumer_resource_data_analysis as cf
import matplotlib.pyplot as plt

folder='structural_stability'
filename='critical_dynamical_syntrophies_NR50_NS25'

alpha_mode=['FC', 'NIS', 'LRI', 'RS']
additional_alpha_mode=['RNISC', 'LNISC', 'NISCC']

alpha_label=['FC', 'NIS', 'LRI', 'RS']+['RNISC', 'LNISC', 'NISCC']
alpha_mode_colours+=['orange', 'pink', 'grey']


file = folder+'/'+filename+'.out'
mat_name=np.loadtxt(file, usecols=0, dtype='U')
data=np.loadtxt(file, usecols=tuple([i for i in range(1,5)]))
# crit_syn dictionary is organized this way
# crit_syn[Nest][Conn][alpha mode]
nest = []
conn = []
for i in range(len(mat_name)):
    a = mat_name[i].split("_",7)
    actual_nest = float(a[5][4:])
    round_nest = cf.closest_element_in_list(actual_nest, all_nestedness)
    actual_conn = float(a[6][4:-4])
    round_conn = cf.closest_element_in_list(actual_conn, all_connectance)
    nest.append(round_nest)
    conn.append(round_conn)
crit_syn = {}

# build structural stability dictionary
for n in all_nestedness:
    dnest = {}
    for c in all_connectance:
        dconn = {}
        index=0
        for al in alpha_mode:
            dalpha={}
            if index < len(data[0]):
                for i in range(len(mat_name)):
                    if nest[i]==n and conn[i]==c:
                        dalpha['value']=data[i,index]
                index+=1
            dconn[al]=dalpha
        dnest[c]=dconn
    crit_syn[n]=dnest

# add additional alpha modes
for al in additional_alpha_mode:
    new_data=folder+'/'+filename+'_'+al+'.out'
    new_matrices=np.loadtxt(new_data, usecols=0, dtype='U')

    nest=[]
    conn=[]

    for i in range(len(new_matrices)):
        a = new_matrices[i].split("_",7)
        actual_nest = float(a[5][4:])
        round_nest = cf.closest_element_in_list(actual_nest, all_nestedness)
        actual_conn = float(a[6][4:-4])
        round_conn = cf.closest_element_in_list(actual_conn, all_connectance)
        nest.append(round_nest)
        conn.append(round_conn)

    local_crit_syn=np.loadtxt(new_data, usecols=1)
    for n in all_nestedness:
        for c in all_connectance:
            crit_syn[n][c][al]={}
            for i in range(len(new_matrices)):
                if nest[i]==n and conn[i]==c:
                    crit_syn[n][c][al]['value']=local_crit_syn[i]

alpha_mode+=additional_alpha_mode

for k in range(len(alpha_mode)):
    alpha = alpha_mode[k]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    title=alpha_label[k]
    ylabel=r'$\alpha_C^D(\gamma_0=0.75, S_0=0.05)$'
    for nest in all_nestedness:
        data = [crit_syn[nest][k][alpha] for k in all_connectance]
        to_plot = np.array([[all_connectance[i],data[i]['value']] for i in range(len(data)) if data[i]])
        ax.plot(to_plot[:,0], to_plot[:,1], label=r'$\eta_G\approx'+str(nest)+'$')
    ax.set_xlabel(r'$\kappa_G$')
    ax.set_ylabel(ylabel)
    ax.legend(bbox_to_anchor=(1.,1.))
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig('plots/critical_dynamical_syntrophy_fixed_nest_'+alpha+'.pdf')
    plt.close()


    fig = plt.figure()
    ax = fig.add_subplot(111)
    for conn in all_connectance:
        data = [crit_syn[eta][conn][alpha] for eta in all_nestedness]
        to_plot = np.array([[all_nestedness[i],data[i]['value']] for i in range(len(data)) if data[i]])
        ax.plot(to_plot[:,0], to_plot[:,1], label=r'$\kappa_G\approx'+str(conn)+'$')
    ax.set_xlabel(r'$\eta_G$')
    ax.set_ylabel(ylabel)
    ax.legend(bbox_to_anchor=(1.,1.))
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig('plots/critical_dynamical_syntrophy_fixed_conn_'+alpha+'.pdf')
    plt.close()

for conn in all_connectance:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ylabel=r'$\alpha_C^D(\gamma_0=0.75, S_0=0.05)$'
    title = r'$\kappa_G='+str(conn)+'$'
    for k in range(len(alpha_mode)):
        data = [crit_syn[eta][conn][alpha_mode[k]] for eta in all_nestedness]
        to_plot = np.array([[all_nestedness[i],data[i]['value']] for i in range(len(data)) if data[i]])
        ax.plot(to_plot[:,0], to_plot[:,1], label=alpha_label[k], color=alpha_mode_colours[k])
    ax.set_xlabel(r'$\eta_G$')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(bbox_to_anchor=(1.,1.))
    fig.tight_layout()
    fig.savefig('plots/critical_dynamical_syntrophy_Conn'+str(conn)+'.pdf')
    plt.close()
