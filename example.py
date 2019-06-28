#!/usr/local/bin/python3
from microbial_community.microbial_model import Model, plot_time_evolution
import matplotlib.pyplot as plt
import numpy as np

# set printing options
np.set_printoptions(precision=2)

# this example creates a random biological system 
# we consider a nested ecological network consisted of "mainly specialists" 
# (i.e. it is nested but species eat mostly from one resource)
NR_ = 3
NS_ = 3

# the following parameters specify how you want to construct your ecosystem

# this is how you choose the parameters of your system (if you don't know what to put, leave it like that)
construction_mode = 'gamma_construction'
# the type of network you want
network_type = 'mainly specialist nested'
# the energy dissipation of each species
dissipation = 'random'

mode = (construction_mode, network_type, dissipation)

# creation of the model
model = Model()

# we generate random parameters (according to the mode we specified), by default we take an equilibrium of the system
print('We generate a model with random parameters : ')
#model.generate_random_parameters(NR_, NS_, mode)


# Now that the model is generated, we can do a couple of things with it
# save it to an external file (no extension)
model.save('a_cool_model')

# load a new one from an external file (no extension either)
model2 = Model()
model2.load('another_cool_model')
print('We load a model with the parameters : ')
model2.print_parameters()


fig1 = plt.figure(figsize=(10,5))
gs = fig1.add_gridspec(1,2)
axis1 = fig1.add_subplot(gs[0,0])
axis2 = fig1.add_subplot(gs[0,1])

fig2 = plt.figure()
gs = fig2.add_gridspec(1,2)
axis3 = fig2.add_subplot(gs[0,0])
axis4 = fig2.add_subplot(gs[0,1])



# we can make it time evolve from a set of initial conditions
R0 = [1., 1., 1.]
S0 = [1., 1., 1.]
sol = model2.time_evolution(R0, S0, np.linspace(0, 100, 1000))
# and plot it on an axis
plot_time_evolution(axis1, sol)

# we can also plot the eigenspectrum of the jacobian 
model2.plot_eigenspectrum(axis2)

# print the jacobian with the current values
Jeq = model2.jacobian()
print('jacobian at equilibrium = ')
print(Jeq)

# or even draw the two food networks associated to the system
model2.draw_food_network(axis3, 'consumption')
model2.draw_food_network(axis4, 'byproduct release')

plt.tight_layout()
plt.show()

