#============ VARIABLES TO GENERATE FIGURE 2 =======================#
plotting_properties=dict({
    'intrashift': 0.1,  # distance between two boxes at the same alpha0
    'intershift': 0.5,  # distance between two boxes (last at alpha0 and first at the next value of alpha0)
    'width': 1,         # box width
    'data file': 'results/all_data_NR25_NS25.out', # location of the data file produced with compute_feasibility_data.py
    'save folder': 'main_figures',
    'mode': 'test'      # test mode for testing stuff out, plot mode for plotting data
})

#============ DO NOT MODIFY BELOW ==================================#
import consumer_resource_data_analysis.data_loading as dl
import consumer_resource_data_analysis.data_plotting as dp
import matplotlib.pyplot as plt
import numpy as np


plotting_properties['alpha mode']=dp.alpha_mode
plotting_properties['alpha0']=dp.alpha0

data_frame = dl.load_data_frame(plotting_properties['data file'])
print(data_frame.columns)
to_plot=['G connectance', 'G nestedness','A connectance','A nestedness', 'A-mode','alpha0','feasible volume']
print(data_frame[data_frame['alpha0']==0][to_plot].sort_values(by=['G connectance', 'G nestedness','A-mode']).to_string())

if plotting_properties['mode']=='plot':
    figA = plt.figure("Figure 2A")
    axA = figA.add_subplot(111)
    axA = dp.plot_figure_2A(axA, data_frame, plotting_properties)
    figA.tight_layout()
    figA.savefig(plotting_properties['save folder']+"/figure_2A.pdf")

    figB = plt.figure("Figure 2B")
    axB = figB.add_subplot(111)
    axB = dp.plot_figure_2B(axB, data_frame, plotting_properties)
    figB.tight_layout()
    figB.savefig(plotting_properties['save folder']+"/figure_2B.pdf")


    figC = plt.figure("Figure 2C")
    axC = figC.add_subplot(111)
    axC = dp.plot_figure_2C(axC, data_frame, plotting_properties)
    figC.tight_layout()
    figC.savefig(plotting_properties['save folder']+"/figure_2C.pdf")


    ### ATTENTION : MODIFIE LE DATA FRAME ###
    # figD = plt.figure("Figure 2D")
    # axD = []
    # string_axs = int(str(len(plotting_properties['alpha mode']))+"2")
    # for i in range(0,len(plotting_properties['alpha mode'])):
    #     int_subplotA = int(str(string_axs)+str(i*2))+1
    #     int_subplotB = int(str(string_axs)+str(i*2))+2
    #     axD.append([figD.add_subplot(int_subplotA), figD.add_subplot(int_subplotB)])
    # axD = dp.plot_against_matrix_properties(axD, data_frame, plotting_properties,'av. dominant eigenvalue')
    # figD.tight_layout()

    plotting_properties['alpha mode']=['random_structure', 'optimized_matrix']
    plotting_properties['alpha0']=np.sort(data_frame['alpha0'].unique())

    figE = plt.figure("Figure 2E")
    axE = figE.add_subplot(111)
    axE = dp.plot_figure_2E(axE, data_frame, plotting_properties)
    figE.tight_layout()
    figE.savefig(plotting_properties['save folder']+"/figure_2E.pdf")

    figF = plt.figure("Figure 2F")
    axF = figF.add_subplot(111)
    axF = dp.plot_figure_2F(axF, data_frame, plotting_properties)
    figF.tight_layout()
    figF.savefig(plotting_properties['save folder']+"/figure_2F.pdf")

    figG = plt.figure("Figure 2G")
    axG = figG.add_subplot(111)
    axG = dp.plot_figure_2G(axG, data_frame, plotting_properties)
    figG.tight_layout()
    figG.savefig(plotting_properties['save folder']+"/figure_2G.pdf")


    plt.show()
