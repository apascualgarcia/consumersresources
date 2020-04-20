import numpy as np
import matplotlib.pyplot as plt

folder='data_output'
file2 = 'eigenvalues_Nest0.4_Conn0.3312_alpha0.0013.out'
file1 = 'eigenvalues_Nest0.4_Conn0.3712_alpha0.0039.out'
data1 = np.loadtxt(folder+'/'+file1)
data2 = np.loadtxt(folder+'/'+file2)
label1 = r'$\kappa=0.3312$'
label2 = r'$\kappa=0.3712$'

pos_data1 = data1[data1 > 0]
negative_data1 = data1[data1 < 0]
zero_data1 = data1[data1==0.]

hist, bins, _ = plt.hist(negative_data1, bins=200)
neg_logbins1 = -np.logspace(np.log10(abs(bins[0])),np.log10(abs(bins[-1])),len(bins))
hist, bins, _ = plt.hist(pos_data1, bins=200)
pos_logbins1 = np.logspace(np.log10(abs(bins[0])),np.log10(abs(bins[-1])),len(bins))


pos_data2 = data2[data2 > 0]
negative_data2 = data2[data2 < 0]
zero_data2 = data1[data2==0.]
hist, bins, _ = plt.hist(negative_data2, bins=200)
neg_logbins2 = -np.logspace(np.log10(abs(bins[0])),np.log10(abs(bins[-1])),len(bins))
plt.close()

fig = plt.figure()
ax = fig.add_subplot(111)
weights1=np.ones_like(data1)/len(data1)
weights2=np.ones_like(data2)/len(data2)

ax.hist(data1, bins=neg_logbins1, weights=weights1, label=label1)
ax.hist(data2, bins=neg_logbins2, weights=weights2, label=label2)
ax.set_xscale('symlog', linthreshx=min(min(abs(neg_logbins1)), min(pos_logbins1)))
ax.set_xlabel(r'$|$Re$(\lambda_1)|$')
ax.legend()
fig.tight_layout()
plt.show()
