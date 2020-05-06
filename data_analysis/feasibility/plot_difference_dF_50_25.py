import numpy as np
import matplotlib.pyplot as plt
import consumer_resource_data_analysis as cf
from scipy.stats import linregress
from consumer_resource_data_analysis import alpha_mode, all_nestedness, all_connectance, label, nest_colours, conn_colours

type='feasibility'
dF_50_filename='feasibility/feasibility_NR50_NS25_full_rank_opt_consumption_mat_NR50_NS25'
dF_25_filename='feasibility/all_mat_feasibility_NR25_NS25_100_points_full_rank_opt_consumption_mat_NR25_NS25'
ylabel=r'$d_F(N_R=25)/d_F(N_R=50)$'

type='local_dynamical_stability'
dF_50_filename='local_dynamical_stability/local_dynamical_stability_NR50_NS25_full_rank_opt_consumption_mat_NR50_NS25'
dF_25_filename='local_dynamical_stability/all_mat_local_dynamical_stability_NR25_NS25_100_points_full_rank_opt_consumption_mat_NR25_NS25'
ylabel=r'$d_D(N_R=25)/d_D(N_R=50)$'


alpha_mode=['fully_connected', 'no_release_when_eat', 'optimal_matrix']

df_50 = []
df_25 = []

for a in alpha_mode:
    name_50 = dF_50_filename+'_decay_rate_'+a+'.out'
    name_25 = dF_25_filename+'_decay_rate_'+a+'.out'

    df_50.append(np.loadtxt(name_50))
    df_25.append(np.loadtxt(name_25))

df_50 = np.array(df_50)
df_25 = np.array(df_25)

# sort data first
# data is organized this way: decay_rate[NR][alpha_mode][nest][conn]
decay_rate={}
dNR25={}
dNR50={}
for i in range(len(alpha_mode)):
    dNR25_almode={}
    dNR50_almode={}

    data_50 = df_50[i]
    data_25 = df_25[i]

    # we first take all NR=50 points and match them to a NR=25 point
    NR = data_25[:,0]
    NS = data_25[:,1]
    nestedness_25=[cf.closest_element_in_list(n, all_nestedness) for n in data_25[:,2]]
    connectance_25=[cf.closest_element_in_list(k, all_connectance) for k in data_25[:,3]]
    decay_rate_25=data_25[:,4]

    nestedness_50=[cf.closest_element_in_list(n, all_nestedness) for n in data_50[:,2]]
    connectance_50=[cf.closest_element_in_list(k, all_connectance) for k in data_50[:,3]]
    decay_rate_50=data_50[:,4]

    # first for NR=25
    for nest in all_nestedness:
        dnest={}
        for conn in all_connectance:
            index_decay_rate=[i for i in range(len(nestedness_25)) if nestedness_25[i]==nest and connectance_25[i]==conn]
            if len(index_decay_rate)>0:
                if len(index_decay_rate)==1:
                    drate=decay_rate_25[index_decay_rate[0]]
                    dnest[conn]=drate
                else:
                    print("Problem, two matrices with the same properties: nest", nest,"conn", conn, "(", alpha_mode[i], ")")
        dNR25_almode[nest]=dnest

    # then for NR=50
    # first for NR=25
    for nest in all_nestedness:
        dnest={}
        for conn in all_connectance:
            index_decay_rate=[i for i in range(len(nestedness_50)) if nestedness_50[i]==nest and connectance_50[i]==conn]
            if len(index_decay_rate)>0:
                if len(index_decay_rate)==1:
                    drate=decay_rate_50[index_decay_rate[0]]
                    dnest[conn]=drate
                else:
                    print("Problem, two matrices with the same properties")
        dNR50_almode[nest]=dnest

    dNR25[alpha_mode[i]]=dNR25_almode
    dNR50[alpha_mode[i]]=dNR50_almode

decay_rate[25]=dNR25
decay_rate[50]=dNR50

# plot df(50) vs df(25)
# for i in range(len(alpha_mode)):
#     points_to_plot=[]
#
#     data_50 = df_50[i]
#     data_25 = df_25[i]
#
#     # we first take all NR=25 points and match them to a NR=50 point
#     NR = data_25[:,0]
#     NS = data_25[:,1]
#     nestedness_25=data_25[:,2]
#     connectance_25=data_25[:,3]
#
#     decay_rate_25=data_25[:,4]
#
#     nestedness_50=data_50[:,2]
#     connectance_50=data_50[:,3]
#
#
#     decay_rate_50=data_50[:,4]
#
#     for k in range(len(connectance_25)) :
#         conn=cf.closest_element_in_list(connectance_25[k], all_connectance)
#         indeces_c=[j for j in range(len(connectance_50)) if cf.closest_element_in_list(connectance_50[j], all_connectance)==conn]
#         nest=cf.closest_element_in_list(nestedness_25[k], all_nestedness)
#         indeces_n=[j for j in range(len(nestedness_50)) if cf.closest_element_in_list(nestedness_50[j], all_nestedness)==nest]
#         common_indices=[j for j in indeces_n if j in indeces_c]
#         if len(common_indices)>0:
#             index=np.min(common_indices)
#             points_to_plot.append([conn, nest, decay_rate_25[k], decay_rate_50[index]])
#     points_to_plot=np.array(points_to_plot)
#     power_law_fit, popt, perr = cf.fit_data(cf.power_function, points_to_plot[:,2], points_to_plot[:,3])
#     slope, intercept, rval, pval, stderr = linregress(points_to_plot[:,2:4])
#
#     # plot the points
#     fig1 = plt.figure()
#     ax = fig1.add_subplot(111)
#     points_identity = np.linspace(min(points_to_plot[:,2]), max(points_to_plot[:,2]), 100)
#     linear_fit = [p*slope+intercept for p in points_identity]
#     power_law_fit = [cf.power_function(x, *popt) for x in points_identity]
#     print(popt)
#     # plot identity as comparison
#     #ax.plot(points_identity, points_identity, marker='None', linestyle='solid')
#     ax.plot(points_identity, power_law_fit, marker='None', linestyle='solid', color='black')
#     ax.set_xlabel(r'$d_F(N_R=25)$')
#     ax.set_ylabel(r'$d_F(N_R=50)$')
#     ax.set_title(label[i])
#     im=ax.scatter(points_to_plot[:,2], points_to_plot[:,3], c=points_to_plot[:,1], cmap='jet', s=50)
#     fig1.colorbar(im, ax=ax, orientation='vertical', label=r'$\eta_G$')
#     fig1.tight_layout()
#     fig1.savefig('plots/decay_rate_25_vs_50_nestedness_'+alpha_mode[i]+'.pdf')
#     plt.close()
#
#     # plot the points
#     fig2 = plt.figure()
#     ax = fig2.add_subplot(111)
#     points_identity = np.linspace(min(points_to_plot[:,2]), max(points_to_plot[:,2]), 100)
#     # plot identity as comparison
#     #ax.plot(points_identity, points_identity, marker='None', linestyle='solid')
#     ax.plot(points_identity, power_law_fit, marker='None', linestyle='solid', color='black')
#     ax.set_xlabel(r'$d_F(N_R=25)$')
#     ax.set_ylabel(r'$d_F(N_R=50)$')
#     ax.set_title(label[i])
#     im=ax.scatter(points_to_plot[:,2], points_to_plot[:,3], c=points_to_plot[:,0], cmap='jet', s=50)
#     fig2.colorbar(im, ax=ax, orientation='vertical', label=r'$\kappa_G$')
#     fig2.tight_layout()
#     fig2.savefig('plots/decay_rate_25_vs_50_connectance_'+alpha_mode[i]+'.pdf')
#     plt.close()
#
#
# plot df(50)/df(25) as a function of the consumption matrix characteristics

for i in range(len(alpha_mode)):
    fig=plt.figure()
    ax=fig.add_subplot(111)
    av_val=[]
    min_conn=0.5
    max_conn=0.4
    for j in range(len(all_nestedness)):
        nest=all_nestedness[j]
        d25=decay_rate[25][alpha_mode[i]][nest]
        d50=decay_rate[50][alpha_mode[i]][nest]
        conn=np.array([k for k in all_connectance if (k in d50) and (k in d25)])
        if np.min(conn) < min_conn:
            min_conn=np.min(conn)
        if np.max(conn) > max_conn:
            max_conn=np.max(conn)
        if(len(conn)>0):
            ratio=[d25[k]/d50[k] for k in conn]
            av_val.append(np.mean(ratio))
            ax.plot(conn, ratio, label=r'$\eta_G\approx'+str(nest)+'$', color=nest_colours[j])
    ax.plot(np.linspace(min_conn, max_conn, 10), [np.mean(av_val) for i in range(10)], marker='', linestyle='solid', alpha=0.7, label='Av. val.', color='black')
    ax.legend(bbox_to_anchor=(1., 1.))
    ax.set_xlabel(r'Connectance $\kappa_G$')
    ax.set_ylabel(ylabel)
    ax.set_title(label[i])
    fig.tight_layout()
    fig.savefig('plots/'+type+'_decay_rate_50_vs_25_fixed_nestedness_'+alpha_mode[i]+'.pdf')

    fig=plt.figure()
    ax=fig.add_subplot(111)
    av_val=[]
    min_nest=0.5
    max_nest=0.4
    for j in range(len(all_connectance)):
        conn=all_connectance[j]
        nestedness=[eta for eta in all_nestedness if (conn in decay_rate[25][alpha_mode[i]][eta]) and (conn in decay_rate[50][alpha_mode[i]][eta])]
        if np.min(nestedness) < min_nest:
            min_nest=np.min(nestedness)
        if np.max(nestedness) > max_nest:
            max_nest=np.max(nestedness)
        ratio=[decay_rate[25][alpha_mode[i]][nest][conn]/decay_rate[50][alpha_mode[i]][nest][conn] for nest in nestedness]
        av_val.append(np.mean(ratio))
        ax.plot(nestedness, ratio, label=r'$\kappa_G\approx'+str(conn)+'$', color=conn_colours[j])
    ax.plot(np.linspace(min_nest, max_nest, 10), [np.mean(av_val) for i in range(10)], marker='', linestyle='solid', alpha=0.7, label='Av. val.', color='black')
    ax.legend(bbox_to_anchor=(1., 1.))
    ax.set_xlabel(r'Ecological overlap $\eta_G$')
    ax.set_ylabel(ylabel)
    ax.set_title(label[i])
    fig.tight_layout()
    fig.savefig('plots/'+type+'_decay_rate_50_vs_25_fixed_connectance_'+alpha_mode[i]+'.pdf')
