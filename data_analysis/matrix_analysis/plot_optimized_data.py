import numpy as np
import pandas as pd
import consumer_resource_data_analysis as cf
import matplotlib.pyplot as plt

##### NON CUSTOMIZABLE PART #######
def plot_optimized_mat_data(data, fig_to_plots):
    if 'nestG-nestA' in fig_to_plots:
        fig1 = plt.figure("Nestedness")
        ax1 = fig1.add_subplot(111, aspect=aspect)
        xlim = [np.min(data['nestG'])-0.1*(np.max(data['nestG'])-np.min(data['nestG'])), 0.1*(np.max(data['nestG'])-np.min(data['nestG']))+np.max(data['nestG'])]
        for i in range(len(cf.all_connectance)):
            connG = cf.all_connectance[i]
            to_plot = data[data['connG']==connG]
            to_plot = to_plot.sort_values('nestG')
            to_plot.plot(x = 'nestG', y = 'nestA', label=r'$\kappa_G='+str(connG)+'$', ax=ax1, xlim = xlim, color=cf.conn_colours[i])
        ax1.set_xlabel(cf.nestedness_label)
        ax1.set_ylabel(r'$\langle \eta_A \rangle $')
        ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        fig1.tight_layout()
        fig1.savefig(save_folder+"optimized_matrices_nestedness_"+save_suffix+'.png', dpi=200)

    if 'connG-connA' in fig_to_plots:
        fig2 = plt.figure("Connectance")
        ax2 = fig2.add_subplot(111, aspect=aspect)
        xlim = [np.min(data['connG'])-0.1*(np.max(data['connG'])-np.min(data['connG'])), 0.1*(np.max(data['connG'])-np.min(data['connG']))+np.max(data['connG'])]
        for i in range(len(cf.all_nestedness)):
            nestG=cf.all_nestedness[i]
            to_plot = data[data['nestG']==nestG]
            to_plot=to_plot.sort_values(by='connG')
            to_plot.plot(x='connG', y='connA', label=r'$\eta_G='+str(nestG)+'$', ax=ax2, xlim=xlim, color=cf.nest_colours[i])
        ax2.set_xlabel(cf.connectance_label)
        ax2.set_ylabel(r'$\langle \kappa_A \rangle$')
        ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        fig2.tight_layout()
        fig2.savefig(save_folder+"optimized_matrices_connectance_"+save_suffix+'.png', dpi=200)

    if 'connG-E' in fig_to_plots:
        fig3 = plt.figure("Energy connectance G")
        ax3 = fig3.add_subplot(111)
        xlim = [np.min(data['connG'])-0.1*(np.max(data['connG'])-np.min(data['connG'])), 0.1*(np.max(data['connG'])-np.min(data['connG']))+np.max(data['connG'])]
        for i in range(len(cf.all_nestedness)):
            nestG = cf.all_nestedness[i]
            to_plot = data[data['nestG']==nestG]
            to_plot = to_plot.sort_values(by='connG')
            to_plot.plot(x='connG', y='E', label=r'$\eta_G='+str(nestG)+'$', ax=ax3, xlim=xlim, color=cf.nest_colours[i])
        ax3.set_xlabel(cf.connectance_label)
        ax3.set_ylabel(r'$\langle E \rangle$')
        ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        fig3.tight_layout()
        fig3.savefig(save_folder+"optimized_matrices_energyvconn_"+save_suffix+'.png', dpi=200)

    if 'nestG-E' in fig_to_plots:
        fig4 = plt.figure("Energy nestedness")
        ax4 = fig4.add_subplot(111)
        xlim = [np.min(data['nestG'])-0.1*(np.max(data['nestG'])-np.min(data['nestG'])), 0.1*(np.max(data['nestG'])-np.min(data['nestG']))+np.max(data['nestG'])]
        for i in range(len(cf.all_connectance)):
            connG = cf.all_connectance[i]
            to_plot = data[data['connG']==connG]
            to_plot = to_plot.sort_values(by='nestG')
            to_plot.plot(x='nestG', y='E', label=r'$\kappa_G='+str(connG)+'$', ax=ax4, xlim=xlim, color=cf.conn_colours[i])
        ax4.set_xlabel(cf.nestedness_label)
        ax4.set_ylabel(r'$\langle E \rangle$')
        ax4.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        fig4.tight_layout()
        fig4.savefig(save_folder+"optimized_matrices_energyvnest_"+save_suffix+'.png', dpi=200)

    if 'connG BM' in fig_to_plots:
        fig = plt.figure("ConnG BM")
        ax = fig.add_subplot(111)
        xlim = [np.min(data['connG'])-0.1*(np.max(data['connG'])-np.min(data['connG'])), 0.1*(np.max(data['connG'])-np.min(data['connG']))+np.max(data['connG'])]
        data = data.sort_values(by='connG')
        data.plot(x='connG', y='nestG', ax=ax, xlim=xlim, label=r'$\langle \eta_G \rangle$')
        data.plot(x='connG', y='connA', ax=ax, xlim=xlim, label=r'$\langle \kappa_A \rangle$')
        data.plot(x='connG', y='nestA', ax=ax, xlim=xlim, label=r'$\langle \eta_A \rangle$')
        ax.set_xlabel(cf.connectance_label)
        ax.set_title('Both G and A modified by the algorithm')
        ax.legend(loc='best')
        fig.tight_layout()
        fig.savefig(save_folder+"optimized_matrices_connGBM_"+save_suffix+'.png', dpi=200)


    return
##### END OF NON CUSTOMIZABLE PART ####


to_plot=[
    # 'connG-connA',
    # 'nestG-nestA',
    # 'connG-E',
    # 'nestG-E',
    'connG BM'
    ]
data_file_path = r"data_output/opt_mat_data_alpha0=1_intra_specific_syntrophy=allowed_both_modified.csv"
#data_file_path = r"data_output/optimized_matrix_data_alpha0=1_intra_specific_syntrophy=allowed_no_iss.csv"
#data_file_path = r"data_output/optimized_matrix_data_alpha0=1_intra_specific_syntrophy=allowed_5millionssteps.csv"
save_suffix = "alpha0=1_bm"
save_folder = "plots/"
# set aspect ratio for connectance and nestedness graphs
aspect = 'auto'
plot_optimized_mat_data(data=pd.read_csv(data_file_path,sep="\t"), fig_to_plots=to_plot)
