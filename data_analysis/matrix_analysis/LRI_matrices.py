import numpy as np
import matplotlib.pyplot as plt

nestedness_gamma=[0.158059,0.104794,0.245902,0.252632,0.20227,0.347044,0.346769,0.354089,0.29744,0.302762,0.446518,0.445459,0.452353,0.454738,0.39726,0.396236,0.404945,0.547425,0.547794,0.546572,0.498845,0.496306,0.504988,0.599869,0.595086,0.596491]
connectance_gamma=[0.1808,0.1296,0.184,0.2336,0.168,0.2208,0.272,0.3216,0.2272,0.2816,0.2672,0.3168,0.384,0.424,0.272,0.3312,0.3712,0.3344,0.3824,0.4176,0.32,0.3824,0.416,0.3168,0.3712,0.4176]

nestedness_alpha_fully_connected=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
nestedness_alpha_no_intraspecific_syntrophy=[0.945043,0.957965,0.945856,0.903195,0.925289,0.940053,0.91143,0.884183,0.908385,0.927579,0.933344,0.894015,0.83829,0.813648,0.929228,0.906773,0.871203,0.918628,0.885086,0.863092,0.899933,0.834159,0.849403,0.910508,0.876804,0.857403]
nestedness_alpha_optimal_matrix=[0.401709,0.423445,0.894784,0.578876,0.896161,0.757812,0.72815,0.70041,0.851432,0.736762,0.885411,0.773136,0.706411,0.60647,0.859623,0.747008,0.683717,0.819983,0.763068,0.723751,0.811792,0.767153,0.693828,0.854336,0.822538,0.778114]

connectance_alpha_fully_connected=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
connectance_alpha_no_intraspecific_syntrophy=[0.9184,0.9392,0.9136,0.8736,0.9008,0.8944,0.8784,0.8448,0.864,0.8896,0.8848,0.848,0.7936,0.7776,0.8752,0.864,0.8368,0.848,0.808,0.8016,0.84,0.7936,0.7856,0.8256,0.8032,0.792]
connectance_alpha_optimal_matrix=[0.1648,0.1312,0.1744,0.224,0.1872,0.2064,0.2592,0.3456,0.224,0.2848,0.248,0.3216,0.3648,0.4544,0.288,0.3536,0.3712,0.3312,0.3856,0.408,0.312,0.4064,0.4224,0.3296,0.3328,0.4032]

nestedness_alpha=[nestedness_alpha_fully_connected, nestedness_alpha_no_intraspecific_syntrophy,nestedness_alpha_optimal_matrix]
connectance_alpha=[connectance_alpha_fully_connected, connectance_alpha_no_intraspecific_syntrophy, connectance_alpha_optimal_matrix]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(connectance_gamma, nestedness_alpha[1], linestyle='', markersize=8, markeredgewidth=2, marker='^', label='intraspecific syntrophy restricted')
ax.plot(connectance_gamma, nestedness_alpha[2], linestyle='', markersize=8, markeredgewidth=2, marker='o', label='optimal LRI')
ax.legend()
ax.set_xlabel(r'$\kappa_G$')
ax.set_ylabel(r'$\eta_A$')
ax.set_title(r'$N_R=25, N_S=25$')
fig.tight_layout()
fig.savefig("plots/structure_alpha_matrix_NR25_N25_with_connectance.pdf")

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.plot(nestedness_gamma, nestedness_alpha[1], linestyle='', markersize=8, markeredgewidth=2, marker='^', label='intraspecific syntrophy restricted')
ax2.plot(nestedness_gamma, nestedness_alpha[2], linestyle='', markersize=8, markeredgewidth=2, marker='o', label='optimal LRI')
ax2.legend()
ax2.set_title(r'$N_R=25, N_S=25$')
ax2.set_xlabel(r'$\eta_G$')
ax2.set_ylabel(r'$\eta_A$')
fig2.tight_layout()
fig2.savefig("plots/structure_alpha_matrix_NR25_NS25_with_nestedness.pdf")
plt.close()


nestedness_gamma=[0.145882,0.154219,0.0971751,0.246764,0.247778,0.254909,0.195578,0.202117,0.345765,0.345419,0.345238,0.354754,0.296258,0.297425,0.304865,0.445187,0.447074,0.445106,0.445314,0.445637,0.453766,0.397492,0.395674,0.396108,0.395657,0.403933,0.545407,0.545063,0.545253,0.545132,0.545487,0.4968,0.495397,0.496476,0.496786,0.495426,0.595393,0.595861,0.595165,0.597376,0.596029]
connectance_gamma=[0.0848,0.12,0.0832,0.1224,0.1656,0.2224,0.1192,0.168,0.1808,0.2296,0.2792,0.332,0.1792,0.2296,0.2792,0.1704,0.2296,0.2768,0.3328,0.3824,0.4336,0.176,0.22,0.2752,0.32,0.3848,0.2176,0.2776,0.332,0.368,0.4344,0.224,0.2736,0.3288,0.3672,0.4272,0.2344,0.268,0.3344,0.3832,0.4272]

nestedness_alpha_fully_connected=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
nestedness_alpha_no_intraspecific_syntrophy=[0.980687,0.966369,0.981784,0.974462,0.962642,0.924326,0.967477,0.945707,0.955671,0.920776,0.900366,0.868611,0.954252,0.919453,0.901451,0.968028,0.938902,0.923329,0.886939,0.844603,0.823433,0.957674,0.945961,0.907218,0.881248,0.872428,0.944447,0.923615,0.914747,0.89324,0.85735,0.945781,0.924918,0.892028,0.880418,0.833122,0.953639,0.946899,0.911663,0.878035,0.860932]
nestedness_alpha_optimal_matrix=[0.511807,0.943503,0.832301,0.85306,0.948687,0.827428,0.906593,0.879521,0.907766,0.912001,0.921549,0.826494,0.910907,0.952337,0.861053,0.92549,0.937541,0.908804,0.89722,0.890082,0.747951,0.885721,0.901171,0.85508,0.939006,0.768711,0.962301,0.965276,0.91303,0.955467,0.929838,0.904709,0.963369,0.945238,0.963869,0.887956,0.940344,0.937178,0.96689,0.944764,0.935489]

connectance_alpha_fully_connected=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
connectance_alpha_no_intraspecific_syntrophy=[0.9496,0.9376,0.964,0.9408,0.928,0.8856,0.9368,0.9128,0.9104,0.8704,0.8552,0.8264,0.9112,0.884,0.8632,0.92,0.8832,0.8656,0.8304,0.7952,0.7784,0.9184,0.8968,0.852,0.828,0.8304,0.8736,0.8576,0.8392,0.8256,0.784,0.888,0.8592,0.8256,0.812,0.788,0.8816,0.8728,0.8328,0.7952,0.7808]
connectance_alpha_optimal_matrix=[0.0832,0.1288,0.0856,0.1312,0.1728,0.24,0.1136,0.1672,0.1848,0.2128,0.2752,0.3424,0.1688,0.232,0.292,0.184,0.2208,0.292,0.3312,0.3688,0.436,0.1592,0.2432,0.2696,0.3312,0.3984,0.2328,0.2848,0.3312,0.3824,0.4192,0.2288,0.276,0.3328,0.3976,0.4168,0.2488,0.2888,0.3264,0.3856,0.4144]

nestedness_alpha=[nestedness_alpha_fully_connected, nestedness_alpha_no_intraspecific_syntrophy,nestedness_alpha_optimal_matrix]
connectance_alpha=[connectance_alpha_fully_connected, connectance_alpha_no_intraspecific_syntrophy, connectance_alpha_optimal_matrix]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(connectance_gamma, nestedness_alpha[1], linestyle='', markersize=8, markeredgewidth=2, marker='^', label='intraspecific syntrophy restricted')
ax.plot(connectance_gamma, nestedness_alpha[2], linestyle='', markersize=8, markeredgewidth=2, marker='o', label='optimal LRI')
ax.legend()
ax.set_xlabel(r'$\kappa_G$')
ax.set_ylabel(r'$\eta_A$')
ax.set_title(r'$N_R=50, N_S=25$')
fig.tight_layout()
fig.savefig("plots/structure_alpha_matrix_NR50_N25_with_connectance.pdf")

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.plot(nestedness_gamma, nestedness_alpha[1], linestyle='', markersize=8, markeredgewidth=2, marker='^', label='intraspecific syntrophy restricted')
ax2.plot(nestedness_gamma, nestedness_alpha[2], linestyle='', markersize=8, markeredgewidth=2, marker='o', label='optimal LRI')
ax2.legend()
ax2.set_title(r'$N_R=50, N_S=25$')
ax2.set_xlabel(r'$\eta_G$')
ax2.set_ylabel(r'$\eta_A$')
fig2.tight_layout()
fig2.savefig("plots/structure_alpha_matrix_NR50_NS25_with_nestedness.pdf")
