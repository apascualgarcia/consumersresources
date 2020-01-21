import matplotlib.pyplot as plt
import numpy as np

plt.rc('lines', marker=None)
plt.rc('lines', linestyle='solid')

NR = 4
NS = 4
folder_path = 'data_output'
filename = 'test'
save_folder='plots/Typical_time_evolution'
added_text='_low_threshold'

def compute_log_derivative(t, f):
    deriv = []
    time = []
    for i in range(len(f)-1):
        time.append(0.5*(t[i+1]+t[i]))
        deriv.append((f[i+1]-f[i])/(f[i]*(t[i+1]-t[i])))
    return np.array(time), np.array(deriv)

data = np.loadtxt(folder_path + "/" + filename +'.out')


t = data[:, 0]
resources = data[:, 1:NR + 1]
consumers = data[:, NR + 1: NS + NR + 1]

fig1 = plt.figure('Time evolution')

ax1 = fig1.add_subplot(211)
ax1.plot(t, resources, color='green')
ax1.set_ylabel(r'Resources')
ax1.xaxis.set_major_formatter(plt.NullFormatter())

ax2 = fig1.add_subplot(212)
ax2.set_xlabel(r'$t$')
ax2.set_ylabel(r'Consumers')
ax2.plot(t, consumers, color='blue')
fig1.tight_layout()
fig1.savefig(save_folder+'/Typical_time_evolution_resources_species'+added_text)

index=0
max_value=10.
for i in range(NS):
    time_deriv, deriv_cons = compute_log_derivative(t, consumers[:, i])
    if(max(abs(deriv_cons[:-1]))>max_value):
        index=i
        max_value = max(abs(deriv_cons))
time_deriv, deriv_cons = compute_log_derivative(t, consumers[:, index])

fig3 = plt.figure('Maximal absolute time log derivative')
ax3 = fig3.add_subplot(111)
ax3.plot(time_deriv, abs(deriv_cons), color='red')
ax3.set_yscale('log')
ax3.set_xlabel(r'$t$')
ax3.set_ylabel(r'$\|d\ln(N_i)/dt\|$')
ax3.set_ylim(bottom=min(abs(deriv_cons)))
fig3.tight_layout()
fig3.savefig(save_folder+'/Typical_time_evolution_log_derivative'+added_text)
