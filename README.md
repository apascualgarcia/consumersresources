This repository is made to study the dynamics of microbial communities.
A part of the research lead so far on this subject is explained in microbial_community_review.pdf

The folder microbial_community is a Python package developed to run simulations of the model considered here.
other scripts contains some unorganized (and not maintained) scripts that were used for testing and the tex file file of microbial_community_review

example.py is a commented Python script that gives examples of the main functions one can use in microbial_community.

Here is however a list of the user accessible functions of microbial_community and its modules:
  * microbial_model : it contains most of the functions related to the Model class.
    - Model.retrieve_model_parameters() : after a model has been initiated (either manually by setting the data of the model by hand or by another function), allows to retrieve the following parameters of the model: (Req, Seq, alpha, gamma, sigma, lambdaa, delta, mu) (as a tuple).
    - Model.generate_random_parameters(NR, NS, param) : allows, after the creation of a model, to set its parameters randomly, according to the way you want (this still needs some proof checking and working). I will update it soon to make sure the different choices work. For now, keep the parameters found in example.py.
    - Model.set_parameters(set_parameters) : allows to manually set the parameters of the model. Input, a dictionary containing the parameters with their different keys (must at least have Req, Seq, alpha, gamma, sigma, lambda, delta, mu for the model to minimally work).
    - Model.print_parameters() : prints on the console the parameters of the Model.
    - Model.draw_food_network(axis, ntype): draws the food network ntype (string either 'consumption' or 'byproduct release') on axis (an axis of a matplotlib figure).
    - Model.jacobian(): returns the jacobian of the model for the current parameters_draw.
    - Model.time_evolution(Rinit, Sinit, t) : returns the time evolution of the model (R(t), S(t), t) on the interval t, assuming that R(t0) = Rinit and S(t0) = Sinit (t0 being the first point of t).
    - Model.save(path) : saves the model as a .npy file on path (path must end with the wanted name, without an extension).
    - Model.load(path): after initiating a model, allows to load its parameters from a model previously saved with Model.save on path.
    - Model.plot_eigenspectrum(axis) : plots the eigenspectrum of the jacobian for the given parameters.
    - plot_time_evolution(axis, sol) : allows to plot on axis the time evolution of the resources and species of the model (we assume sol is the output of Model.time_evolution).

  * parameters_draw : internal module used to draw the parameters for the creation of the model. It should not be called by the user
