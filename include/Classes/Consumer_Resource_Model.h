#ifndef CONSUMER_RESOURCE_MODEL_H
#define CONSUMER_RESOURCE_MODEL_H

#include "Custom_types.h"
#include "Model_Parameters.h"
#include "Metaparameters.h"
#include "Dynamical_variables.h"
#include "Extinction.h"
#include <iostream>
#include <string>

class CRModel{
private:
  Metaparameters* metaparameters;
  ntensor* eq_vals; // we allow the possibility of multiple equilibria (hence vector of vector)
  Model_parameters* model_param;

public:
  CRModel();
  CRModel(Metaparameters&);
  CRModel(const foodmatrix&, Metaparameters&);
  CRModel(Model_parameters*);
  void attempt_to_build_model(const foodmatrix&,Metaparameters&);
  nvector equations_of_evolution(const Dynamical_variables&) const; // returns the value of the RHS of the equations of evolution
  nmatrix jacobian_at_equilibrium() const;
  ncvector eigenvalues_at_equilibrium() const;
  nmatrix jacobian(const Dynamical_variables&) const; // returns the jacobian for the given dynamical variables
  void save(std::ostream&) const; // outputs the model to the external file
  std::ostream& display(std::ostream&) const;
  bool energy_constraint() const;
  bool constraints_fulfilled(const Metaparameters& m) const;
  bool positive_parameters() const;
  bool dynamically_stable() const;
  void save_simulation() const;
  void save_jacobian_at_equilibrium(std::string) const;
  void write_time_evolution(const Dynamical_variables&, ntype) const;
  void write_time_evolution_from_equilibrium() const;
  void write_time_evolution_until_equilibrium(const Dynamical_variables &, ntype, ntype) const;
  void write_death_rates(std::string) const;
  nmatrix time_evolution(const Dynamical_variables&, ntype) const ;
  Dynamical_variables perturb_equilibrium() const;
  void perturb_parameters() const;
  void perturb_parameters(const ntype &) const;
  void save_new_equilibrium(const Extinction&) const;
  double get_m0() const;
  double get_d0() const;
  nvector get_m() const;
  nvector get_d() const;
  ntype get_resilience_jacobian() const;
  ntype get_resilience_dynamical_stability(const ntype& delta=0.);

  bool has_linearly_stable_eq() const;

  nmatrix get_first_equilibrium() const;

  /* returns dynamical variables perturbed away from their equilibrium value*/
  nmatrix perturb_abundances(const ntype& );

  /* returns the extinction properties with the initial values of abundances */
  Extinction evolve_until_equilibrium_from_abundances(const nmatrix& , ntype threshold=1e-9, eqmode eq_mode = oneextinct) const;
  /* returns the evolution from equilibrium (assumes we are not at equilibrium with R^* and S^*) -> finds the new eq*/
  Extinction evolve_until_equilibrium(ntype, eqmode eq_mode=convergence) const;
  /* returns the general extinction properties for the initial values init_val */
  Extinction evolve_until_equilibrium_general(const nmatrix& init_val, ntype threshold, eqmode eq_mode) const;
};

#endif
