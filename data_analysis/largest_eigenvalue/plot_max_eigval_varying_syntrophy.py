import numpy as np
import matplotlib.pyplot as plt
import consumer_resource_data_analysis as cf
import os
import matplotlib.tri as tr
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit
from matplotlib.colors import LogNorm
from consumer_resource_data_analysis import alpha_mode, alpha_mode_colours,label, alpha0, all_nestedness, all_connectance
import copy

filename = 'largest_eigenvalue/largest_eigenvalue_NR25_NS25_100_points_full_rank_opt_consumption_mat_NR25_NS25'
Npoints = 10000
square_size=8


cmap = plt.cm.get_cmap('jet_r')
colors = [cmap(i/10) for i in range(len(alpha0))]

# largest_eigenvalue region[alpha_mode][alpha0][connectance][nestedness][gamma0][S0] contains the largest_eigenvalue of said point
largest_eigenvalue_region = []
for al_mo in alpha_mode:
    local_vector=[]
    for a in alpha0:
        file = filename+'_'+al_mo+'_optimal_LRI_alpha0='+str(a)
        cf.remove_strings_from_file('optimal_matrices/consumption/Nr25_Nc25', file)
        file = filename+'_'+al_mo+'_optimal_LRI_alpha0='+str(a)+'_filtered.out'
        local_data=np.loadtxt(file, dtype=complex)
        local_vector.append(local_data)
    largest_eigenvalue_region.append(local_vector)

largest_eigenvalue_region = np.array(largest_eigenvalue_region)
connectance = np.real(largest_eigenvalue_region[0,0][:,3])
nestedness = np.real(largest_eigenvalue_region[0,0][:,2])

alpha0=np.array(alpha0)

# plot largest eigenvalue plot at no syntrophy for all matrices
for k in range(len(largest_eigenvalue_region[0,0])):
    for j in range(len(largest_eigenvalue_region[0])):
        data = np.ma.masked_invalid(largest_eigenvalue_region[:,j,k])
        if(not(data[:,6::3].mask.all())):
            fig, axs = plt.subplots(1, len(largest_eigenvalue_region), sharey=True, sharex=True, figsize=(3.5*len(largest_eigenvalue_region),4.5))
            max_ev=-1000
            min_ev=1000
            # first find the largest and smallest eigenvalues for the three regimes
            for i in range(len(largest_eigenvalue_region)):
                largest_ev=data[i,6::3]
                if(np.min(np.real(largest_ev)) < min_ev):
                    min_ev = np.min(np.real(largest_ev))
                if(np.max(np.real(largest_ev)) > max_ev):
                    max_ev = np.max(np.real(largest_ev))

            for i in range(len(largest_eigenvalue_region)):
                gamma0=np.real(data[i,4::3])
                S0=np.real(data[i,5::3])
                NR=np.real(data[i,0])
                NS=np.real(data[i,1])
                nestedness=np.real(data[i,2])
                connectance=np.real(data[i,3])
                largest_ev=data[i,6::3]


                axs[i].set_aspect('equal')
                im=axs[i].scatter(gamma0, S0, c=np.real(largest_ev), s=square_size, marker='s', vmin=min_ev, vmax=max_ev, cmap='jet')
                axs[i].set_xlabel(r'$\gamma_0$')
                axs[i].set_xticks([0, 0.5, 1])
                axs[i].set_xticklabels([0, 0.5, 1])
                axs[i].set_title(label[i])

                axs[i].set_xlim(0.,1.)
                axs[i].set_ylim(0.,1.)
            axs[0].set_ylabel(r'$S_0$')
            axs[0].set_yticks([0, 0.5, 1])
            axs[0].set_yticklabels([0, 0.5, 1])
            save_name='NR'+str(int(NR))+'_NS'+str(int(NS))+'_Nest'+str(nestedness)+'_Conn'+str(connectance)+'_alpha0='+str(alpha0[j])
            fig.subplots_adjust(bottom=0.2, top=0.95)
            fig.savefig('plots/largest_eigenvalue_wt_wc_'+save_name+'.pdf')
            cbar_ax = fig.add_axes([0.125, 0.15, 0.75, 0.02])
            cbar=fig.colorbar(im, cax=cbar_ax, orientation='horizontal', format='%.1e')
            cbar.set_label(r'$\langle$Re$(\lambda_1)\rangle$')
            fig.savefig('plots/largest_eigenvalue_wt_'+save_name+'.pdf')
            title = r'Re$(\lambda_1)$ for $N_R='+str(int(NR))+', N_S='+str(int(NS))+', \kappa_G='+str(round(connectance,2))\
                        +', \eta_G='+str(nestedness)+'$'
            title +=r' at $\alpha_0='+str(alpha0[j])+'$'
            fig.suptitle(title)
            fig.savefig('plots/largest_eigenvalue_'+save_name+'.pdf')
            plt.close()

critical_alpha0 = []
for j in range(len(largest_eigenvalue_region[0,0])):
    decline_volume = []
    indices_max_vol=[]
    # find largest alpha0 that has non-zero volume
    for k in range(len(largest_eigenvalue_region)):
        volumes=[]
        max_vol = 0
        i = 0
        continue_loop=True
        while continue_loop:
            data = np.ma.masked_invalid(largest_eigenvalue_region[k,i,j])
            gamma0=np.real(data[4::3])
            S0=np.real(data[5::3])
            NR=np.real(data[0])
            NS=np.real(data[1])
            nestedness=np.real(data[2])
            connectance=np.real(data[3])
            largest_ev=data[6::3]

            if not(data[6::3].mask.all()):
                max_vol +=1
            else:
                continue_loop=False
            continue_loop= (continue_loop and (i < len(largest_eigenvalue_region[0])-1))
            i+=1


        for i in range(max_vol):
            data = np.ma.masked_invalid(largest_eigenvalue_region[k,i,j])
            volumes.append(data[6::3].count()/Npoints)
        decline_volume.append(volumes)
        indices_max_vol.append(max_vol)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    local_crit=[]
    for k in range(len(largest_eigenvalue_region)):
        ax.plot(alpha0[:indices_max_vol[k]], decline_volume[k], label=label[k],markersize=10, linewidth=2.5, markeredgewidth=3)
        max_vol = indices_max_vol[k]
        start_fit=max_vol-4
        fit_function = cf.linear_function
        fit_volumes, popt, perror = cf.fit_data(fit_function, alpha0[start_fit:max_vol], decline_volume[k][start_fit:max_vol])
        alpha0_crit, alpha0_error = cf.zero_from_fit(fit_function, popt, perror)
        print(alpha0_crit, alpha0_error)
        local_crit.append(alpha0_crit)
        ax.plot(alpha0[start_fit:max_vol], fit_volumes, marker='None', linestyle='solid')
    critical_alpha0.append(local_crit)
    ax.legend()
    save_name ='NR'+str(int(NR))+'_NS'+str(int(NS))+'_Nest'+str(nestedness)+'_Conn'+str(connectance)
    title=r'$N_R='+str(int(NR))+', N_S='+str(int(NS))+', \eta='+str(nestedness)+', \kappa='+str(connectance)+'$'
    ax.set_title(title)
    ax.set_xlabel(r'$\alpha_0$')
    ax.set_ylabel(r'Vol$\left(\mathcal{D}^{G,A}_{L,1}(\alpha_0)\right)$ (normalized)')
    fig.tight_layout()
    fig.savefig("plots/size_of_dynamical_volume_"+save_name+".pdf")
    #plt.show()
    plt.close()
critical_alpha0=np.transpose(critical_alpha0)

connectance = np.real(largest_eigenvalue_region[0,0][:,3])
nestedness = np.real(largest_eigenvalue_region[0,0][:,2])
for k in range(len(largest_eigenvalue_region)):
    # plot critical alpha0
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for nest in all_nestedness:
        indices = [int(i) for i in range(len(nestedness)) if cf.closest_element_in_list(nestedness[i], all_nestedness)==nest]
        sorted_indices=np.array([indices[a] for a in np.argsort(connectance[indices])])
        connectance_to_plot=[connectance[i] for i in sorted_indices]
        critical_alpha0_to_plot=[critical_alpha0[k][i] for i in sorted_indices]
        ax.plot(connectance_to_plot, critical_alpha0_to_plot, label=r'$\eta_G\approx'+str(nest)+'$', markersize=10, linewidth=2.5, markeredgewidth=3)
    ax.set_xlabel(r'Connectance $\kappa_G$')
    ax.set_ylabel(r'$\alpha_0^D(G,A)$')
    ax.set_title(label[k])
    ax.legend(bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout()
    fig.savefig('plots/largest_eigenvalue_NR'+str(int(NR))+'_NS'+str(int(NS))+'_critical_alpha0_fixed_nestedness_'+alpha_mode[k]+'.pdf')
    plt.close()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for conn in all_connectance:
        indices = [int(i) for i in range(len(connectance)) if cf.closest_element_in_list(connectance[i], all_connectance)==conn]
        sorted_indices=np.array([indices[a] for a in np.argsort(nestedness[indices])])
        nestedness_to_plot=[nestedness[i] for i in sorted_indices]
        critical_alpha0_to_plot=[critical_alpha0[k][i] for i in sorted_indices]
        ax.plot(nestedness_to_plot, critical_alpha0_to_plot, label=r'$\kappa_G\approx'+str(conn)+'$',markersize=10, linewidth=2.5, markeredgewidth=3)
    ax.set_xlabel(r'Ecological overlap $\eta_G$')
    ax.set_ylabel(r'$\alpha_0^D(G,A)$')
    ax.set_title(label[k])
    ax.legend(bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout()
    fig.savefig('plots/largest_eigenvalue_NR'+str(int(NR))+'_NS'+str(int(NS))+'_critical_alpha0_fixed_connectance_'+alpha_mode[k]+'.pdf')
    plt.close()

# now plot largest eigenvalue observed
exponents=[]
for k in range(len(largest_eigenvalue_region[0,0])):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    local_exponents=[]
    for i in range(len(largest_eigenvalue_region)):
        data=np.ma.masked_invalid(largest_eigenvalue_region[i,:,k])
        NR=np.real(data[0, 0])
        NS=np.real(data[0, 1])
        nestedness_=np.real(data[0,2])
        connectance_=np.real(data[0,3])
        largest_eigenvalue=np.real(data[:, 6::3])
        largest_ev=np.max(largest_eigenvalue, axis=1)
        indices_fit=~largest_ev.mask
        fit, popt, error = cf.fit_data(cf.power_function, alpha0[indices_fit], np.abs(largest_ev[indices_fit]))
        local_exponents.append(popt[1])
        ax.plot(alpha0, np.abs(largest_ev), label=label[i],markersize=10, linewidth=2.5, markeredgewidth=3)
        #ax.plot(alpha0[indices_fit], fit, marker='', linestyle='solid')
    exponents.append(local_exponents)
    save_name ='NR'+str(int(NR))+'_NS'+str(int(NS))+'_Nest'+str(nestedness_)+'_Conn'+str(connectance_)
    ax.set_xlim(0, alpha0[-1]*(1.01))
    ax.set_xlabel(r'$\alpha_0$')
    ax.set_ylabel(r'$\max_{(\gamma_0, S_0)\in [0,1]^2}|\langle$Re($\lambda_1$)$\rangle|$')
    ax.set_yscale('linear')
    ax.ticklabel_format(axis="both", style="sci", scilimits=(-2,2))
    ax.legend()
    fig.tight_layout()
    fig.savefig('plots/largest_eigenvalue_varying_syntrophy_'+save_name+'.pdf')
#    plt.show()
    plt.close()
exponents=np.transpose(exponents)

# plot critical exponents that tells you how deeper your eigenvalues go
for k in range(len(largest_eigenvalue_region)):
    # plot exponent
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for nest in all_nestedness:
        indices = [int(i) for i in range(len(nestedness)) if cf.closest_element_in_list(nestedness[i], all_nestedness)==nest]
        sorted_indices=np.array([indices[a] for a in np.argsort(connectance[indices])])
        connectance_to_plot=[connectance[i] for i in sorted_indices]
        exp_to_plot=[exponents[k][i] for i in sorted_indices]
        ax.plot(connectance_to_plot, exp_to_plot, label=r'$\eta\approx'+str(nest)+'$', markersize=10, linewidth=2.5, markeredgewidth=3)
    ax.set_xlabel(r'Connectance $\kappa$')
    ax.set_ylabel(r'Fit exponent')
    ax.set_title(label[k])
    ax.legend(bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout()
    fig.savefig('plots/fit_exponent_largest_eigenvalue_NR'+str(int(NR))+'_NS'+str(int(NS))+'_critical_alpha0_fixed_nestedness_'+alpha_mode[k]+'.pdf')
    plt.close()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for conn in all_connectance:
        indices = [int(i) for i in range(len(connectance)) if cf.closest_element_in_list(connectance[i], all_connectance)==conn]
        sorted_indices=np.array([indices[a] for a in np.argsort(nestedness[indices])])
        nestedness_to_plot=[nestedness[i] for i in sorted_indices]
        exp_to_plot=[exponents[k][i] for i in sorted_indices]
        ax.plot(nestedness_to_plot, exp_to_plot, label=r'$\kappa\approx'+str(conn)+'$',markersize=10, linewidth=2.5, markeredgewidth=3)
    ax.set_xlabel(r'Ecological overlap $\eta$')
    ax.set_ylabel(r'Fit exponent')
    ax.set_title(label[k])
    ax.legend(bbox_to_anchor=(1.0, 1.0))
    fig.tight_layout()
    fig.savefig('plots/fit_exponent_largest_eigenvalue_NR'+str(int(NR))+'_NS'+str(int(NS))+'_critical_alpha0_fixed_connectance_'+alpha_mode[k]+'.pdf')
    plt.close()
