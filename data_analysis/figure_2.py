#============ VARIABLES TO GENERATE FIGURE 2 =======================#
plotting_properties=dict({
    'intrashift': 0.1,  # distance between two boxes at the same alpha0
    'intershift': 0.5,  # distance between two boxes (last at alpha0 and first at the next value of alpha0)
    'width': 1,         # box width
    'data file': 'results/all_data_NR25_NS25_verbose-level=1.out', # location of the data file produced with compute_feasibility_data.py
    'save name': 'main_figures/feasibility/feasible_volume_Nr25_Nc25.pdf',
})

#============ DO NOT MODIFY BELOW ==================================#
import consumer_resource_data_analysis.data_loading as dl
import consumer_resource_data_analysis.data_plotting as dp
import matplotlib.pyplot as plt


plotting_properties['alpha mode']=dp.alpha_mode
plotting_properties['alpha0']=dp.alpha0

data_frame = dl.load_data_frame(plotting_properties['data file'])

figA = plt.figure("Figure 2A")
axA = figA.add_subplot(111)
axA = dp.plot_figure_2A(axA, data_frame, plotting_properties)
figA.tight_layout()

figB = plt.figure("Figure 2B")
axB = figB.add_subplot(111)
axB = dp.plot_figure_2B(axB, data_frame, plotting_properties)
figB.tight_layout()

figC = plt.figure("Figure 2C")
axC = figC.add_subplot(111)
axC = dp.plot_figure_2C(axC, data_frame, plotting_properties)
figC.tight_layout()

plt.show()
