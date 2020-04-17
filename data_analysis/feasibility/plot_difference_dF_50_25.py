import numpy as np
import matplotlib.pyplot as plt
import consumer_resource_data_analysis as cf
from scipy.stats import linregress


dF_50_filename='feasibility/feasibility_NR50_NS25_full_rank_opt_consumption_mat_NR50_NS25'
dF_25_filename='feasibility/feasibility_NR25_NS25_full_rank_opt_consumption_mat_NR25_NS25'

all_nestedness=[0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]
all_connectance=[0.08, 0.13, 0.18, 0.23, 0.28, 0.33, 0.38, 0.43]
alpha_mode=['fully_connected', 'no_release_when_eat', 'optimal_matrix']
label=['fully connected', 'no intraspecific syntrophy', 'LRI matrix']

df_50 = []
df_25 = []

for a in alpha_mode:
    name_50 = dF_50_filename+'_decay_rate_'+a+'.out'
    name_25 = dF_25_filename+'_decay_rate_'+a+'.out'

    df_50.append(np.loadtxt(name_50))
    df_25.append(np.loadtxt(name_25))

df_50 = np.array(df_50)
df_25 = np.array(df_25)

for i in range(len(alpha_mode)):
    points_to_plot=[]

    data_50 = df_50[i]
    data_25 = df_25[i]

    # we first take all NR=25 points and match them to a NR=50 point
    NR = data_25[:,0]
    NS = data_25[:,1]
    nestedness_25=data_25[:,2]
    connectance_25=data_25[:,3]

    decay_rate_25=data_25[:,4]

    nestedness_50=data_50[:,2]
    connectance_50=data_50[:,3]


    decay_rate_50=data_50[:,4]

    for k in range(len(connectance_25)) :
        conn=cf.closest_element_in_list(connectance_25[k], all_connectance)
        indeces_c=[j for j in range(len(connectance_50)) if cf.closest_element_in_list(connectance_50[j], all_connectance)==conn]
        nest=cf.closest_element_in_list(nestedness_25[k], all_nestedness)
        indeces_n=[j for j in range(len(nestedness_50)) if cf.closest_element_in_list(nestedness_50[j], all_nestedness)==nest]
        common_indices=[j for j in indeces_n if j in indeces_c]
        if len(common_indices)>0:
            index=np.min(common_indices)
            points_to_plot.append([conn, nest, decay_rate_25[k], decay_rate_50[index]])
    points_to_plot=np.array(points_to_plot)
    power_law_fit, popt, perr = cf.fit_data(cf.power_function, points_to_plot[:,2], points_to_plot[:,3])
    slope, intercept, rval, pval, stderr = linregress(points_to_plot[:,2:4])

    # plot the points
    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    points_identity = np.linspace(min(points_to_plot[:,2]), max(points_to_plot[:,2]), 100)
    linear_fit = [p*slope+intercept for p in points_identity]
    power_law_fit = [cf.power_function(x, *popt) for x in points_identity]
    print(popt)
    # plot identity as comparison
    #ax.plot(points_identity, points_identity, marker='None', linestyle='solid')
    ax.plot(points_identity, power_law_fit, marker='None', linestyle='solid', color='black')
    ax.set_xlabel(r'$d_F(N_R=25)$')
    ax.set_ylabel(r'$d_F(N_R=50)$')
    ax.set_title(label[i])
    im=ax.scatter(points_to_plot[:,2], points_to_plot[:,3], c=points_to_plot[:,1], cmap='jet', s=50)
    fig1.colorbar(im, ax=ax, orientation='vertical', label=r'$\eta_G$')
    fig1.tight_layout()
    fig1.savefig('plots/decay_rate_25_vs_50_nestedness_'+alpha_mode[i]+'.pdf')
    plt.close()

    # plot the points
    fig2 = plt.figure()
    ax = fig2.add_subplot(111)
    points_identity = np.linspace(min(points_to_plot[:,2]), max(points_to_plot[:,2]), 100)
    # plot identity as comparison
    #ax.plot(points_identity, points_identity, marker='None', linestyle='solid')
    ax.plot(points_identity, power_law_fit, marker='None', linestyle='solid', color='black')
    ax.set_xlabel(r'$d_F(N_R=25)$')
    ax.set_ylabel(r'$d_F(N_R=50)$')
    ax.set_title(label[i])
    im=ax.scatter(points_to_plot[:,2], points_to_plot[:,3], c=points_to_plot[:,0], cmap='jet', s=50)
    fig2.colorbar(im, ax=ax, orientation='vertical', label=r'$\kappa_G$')
    fig2.tight_layout()
    fig2.savefig('plots/decay_rate_25_vs_50_connectance_'+alpha_mode[i]+'.pdf')
    plt.close()
