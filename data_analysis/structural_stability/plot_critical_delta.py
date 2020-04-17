import numpy as np
from consumer_resource_data_analysis import alpha_mode, all_nestedness, all_connectance, alpha_mode_colours
import consumer_resource_data_analysis as cf
import matplotlib.pyplot as plt
from consumer_resource_data_analysis import label as alpha_label

folder='structural_stability'
filename='data_all_structural_stability'
syntrophy_mode=['common_max_syntrophy','own_max_syntrophy', 'no_syntrophy']
syntrophy_label=['Common maximum syntrophy', 'Own individual syntrophy', 'No syntrophy']
syntrophy_marker=['^', 'o', 's']

file = folder+'/'+filename+'.out'
mat_name=np.loadtxt(file, usecols=0, dtype='U')
data=np.loadtxt(file, usecols=tuple([i for i in range(1,19)]))
# str_stab dictionary is organized this way
# str_stab[Nest][Conn][syntrophy mode][alpha mode]
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
str_stab = {}
# build structural stability dictionary
for n in all_nestedness:
    dnest = {}
    for c in all_connectance:
        dconn = {}
        index=0
        for syn in syntrophy_mode:
            dsyn={}
            for al in alpha_mode:
                dalpha={}
                if index < len(data[0]):
                    for i in range(len(mat_name)):
                        if nest[i]==n and conn[i]==c:
                            dalpha['value']=data[i,index]
                            dalpha['error']=data[i,index+1]
                    index+=2
                    dsyn[al]=dalpha
                dconn[syn]=dsyn
        dnest[c]=dconn
    str_stab[n]=dnest
fig_index=0
for m in range(len(syntrophy_mode)):
    for n in range(len(alpha_mode)):
        syn = syntrophy_mode[m]
        alpha = alpha_mode[n]
        if not(syn=='no_syntrophy' and alpha!='fully_connected'):
            title = syntrophy_label[m]+', '+alpha_label[n]
            save_string=syn+'_'+alpha
            fig = plt.figure(fig_index)
            ax = fig.add_subplot(111)
            for nest in all_nestedness:
                data = [str_stab[nest][k][syn][alpha] for k in all_connectance]
                to_plot = np.array([[all_connectance[i],data[i]['value'], data[i]['error']] for i in range(len(data)) if data[i]])
                ax.plot(to_plot[:,0], to_plot[:,1], label=r'$\eta_G\approx'+str(nest)+'$')
            ax.set_ylabel(r'$\Delta_S^*$')
            ax.set_xlabel(r'$\kappa_G$')
            ax.set_title(title)
            ax.legend(bbox_to_anchor=(1.0, 1.0))
            fig.tight_layout()
            fig.savefig('plots/critical_delta_str_stab_fixed_nest_'+save_string+'.pdf')
            fig_index+=1
            plt.close()

            fig = plt.figure(fig_index)
            ax = fig.add_subplot(111)
            for conn in all_connectance:
                data = [str_stab[eta][conn][syn][alpha] for eta in all_nestedness]
                to_plot = np.array([[all_nestedness[i],data[i]['value'], data[i]['error']] for i in range(len(data)) if data[i]])
                ax.plot(to_plot[:,0], to_plot[:,1], label=r'$\kappa_G\approx'+str(conn)+'$')
            ax.set_ylabel(r'$\Delta_S^*$')
            ax.set_xlabel(r'$\eta_G$')
            ax.set_title(title)
            ax.legend(bbox_to_anchor=(1.0, 1.0))
            fig.tight_layout()
            fig.savefig('plots/critical_delta_str_stab_fixed_conn_'+save_string+'.pdf')
            fig_index+=1
            plt.close()
for conn in all_connectance:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    null=[str_stab[eta][conn]['no_syntrophy']['fully_connected'] for eta in all_nestedness]
    title=r'$\kappa_G='+str(conn)+'$'
    for m in range(len(syntrophy_mode)):
        for n in range(len(alpha_mode)):
            syn = syntrophy_mode[m]
            alpha = alpha_mode[n]
            if not(syn=='no_syntrophy' and alpha!='fully_connected'):
                label = syntrophy_label[m]+', '+alpha_label[n]
                data=[str_stab[eta][conn][syn][alpha] for eta in all_nestedness]
                to_plot = np.array([[all_nestedness[i],data[i]['value']-null[i]['value'], data[i]['error']+null[i]['value']] for i in range(len(data)) if data[i]])
                ax.plot(to_plot[:,0], to_plot[:,1], label=label, color=alpha_mode_colours[n], marker=syntrophy_marker[m], markersize=6)
    ax.set_ylabel(r'$\Delta_S(m,G,A)-\Delta_S(m, G, 0)$')
    ax.set_xlabel(r'$\eta_G$')
    ax.set_title(title)
    fig.tight_layout()
    plt.show()
