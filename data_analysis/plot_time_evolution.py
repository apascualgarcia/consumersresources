import matplotlib.pyplot as plt
import numpy as np


NR = 25
NS = 25
folder_path = "/Users/Shared/Cloud/Master/Master_Thesis/Code/consumersresources/TEST_FOLDER/DataOutput"
filename = "test.out"

data = np.loadtxt(folder_path + "/" + filename)


t = data[:, 0]
resources = data[:, 1:NR + 1]
consumers = data[:, NR + 1: NS + NR + 1]
conv_coeff = data[:, -1]

print(consumers)
print(resources)

fig1 = plt.figure('Time evolution resources')
ax1 = fig1.add_subplot(111)
ax1.plot(t, resources, color='green')

fig2 = plt.figure('Time evolution consumers')
ax2 = fig2.add_subplot(111)
ax2.plot(t, consumers, color='blue')

fig3 = plt.figure('Time evolution convergence')
ax3 = fig3.add_subplot(111)
ax3.plot(t, conv_coeff, color='red')
ax3.set_yscale('log')

plt.show()
plt.clf()
