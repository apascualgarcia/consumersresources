import matplotlib.pyplot as plt
import numpy as np

plt.rc('lines', marker=None)
plt.rc('lines', linestyle='solid')

NR = 25
NS = 25
folder_path = 'data_output'
file_ext='.out'
filename = 'test'
save_folder='plots/Problematic_time_evolution'
added_text=''
INTEGRATOR_ZERO=1e-14
title= r"$\Delta=0.01, S_0=0.1, \gamma_0=0.1$ ($\eta=0.35, \kappa=0.08$, av. conv.time=7.62s)"
save_name="8_slow_convergence"


def compute_log_derivative(t, f):
    deriv = []
    time = []
    for i in range(len(f)-1):
        time.append(0.5*(t[i+1]+t[i]))
        deriv.append((f[i+1]-f[i])/(f[i]*(t[i+1]-t[i])))
    return np.array(time), np.array(deriv)

data = np.loadtxt(folder_path + "/" + filename +file_ext)

t = data[:, 0]
steps = range(1,len(t)+1)
resources = data[:, 1:NR + 1]
consumers = data[:, NR + 1: NS + NR + 1]
error = data[:, -1]

fig1 = plt.figure('Timestep evolution')
ax11 = fig1.add_subplot(211)
ax11.plot(steps, resources, color='green')
ax11.set_ylabel(r'Resources')
ax11.set_yscale('linear')
ax11.set_title(title, fontsize=8)
ax11.xaxis.set_major_formatter(plt.NullFormatter())

ax12 = fig1.add_subplot(212)
ax12.set_xlabel(r'Timestep')
ax12.set_ylabel(r'Consumers')
ax12.plot(steps, consumers, color='blue')
ax12.set_yscale('linear')
fig1.tight_layout()
fig1.savefig(save_folder+'/'+save_name+'_resources_species_timestep'+added_text)

fig2 = plt.figure('Timestep evolution log')
ax21 = fig2.add_subplot(111)
ax21.plot(steps, consumers, color='blue')
ax21.set_xlabel(r'Timesteps')
ax21.set_ylabel(r'Consumers')
ax21.set_yscale('log')
ax21.set_title(title,fontsize=8)
fig2.tight_layout()
fig2.savefig(save_folder+'/'+save_name+'_species_log'+added_text)


fig3 = plt.figure('Time evolution')
ax31 = fig3.add_subplot(211)
ax31.plot(t, resources, color='green')
ax31.set_title(title,fontsize=8)
ax31.set_ylabel(r'Resources')
ax31.xaxis.set_major_formatter(plt.NullFormatter())

ax32 = fig3.add_subplot(212)
ax32.set_xlabel(r'Real time')
ax32.set_ylabel(r'Consumers')
ax32.plot(t, consumers, color='blue')
fig3.tight_layout()
fig3.savefig(save_folder+'/'+save_name+'_resources_species_real_time'+added_text)

# index=0
# max_value=10.
# for i in range(NS):
#     if(min(consumers[:, i]) > INTEGRATOR_ZERO):
#         time_deriv, deriv_cons = compute_log_derivative(t, consumers[:, i])
#         if(max(abs(deriv_cons[:-1]))>max_value):
#             index=i
#             max_value = max(abs(deriv_cons))
# time_deriv, deriv_cons = compute_log_derivative(t, consumers[:, index])
# steps_time_deriv = range(1, len(time_deriv)+1)
#
# fig4 = plt.figure('Maximal absolute time log derivative')
# ax41 = fig4.add_subplot(111)
# ax41.set_title(title)
# ax41.plot(steps_time_deriv, abs(deriv_cons), color='red')
# ax41.set_yscale('log')
# ax41.set_xlabel(r'Timestep')
# ax41.set_ylabel(r'$\|d\ln(N_i)/dt\|$')
# ax41.set_ylim(bottom=min(abs(deriv_cons)))
# fig4.tight_layout()
# fig4.savefig(save_folder+'/Typical_time_evolution_log_derivative'+added_text)
