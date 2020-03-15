#ifndef CONSUMER_RESOURCE_MODEL_H
#define CONSUMER_RESOURCE_MODEL_H

#include "Custom_types.h"
#include "Model_Parameters.h"
#include "Metaparameters.h"
#include "Dynamical_variables.h"
#include "Extinction.h"
#include "../Functions/model_characteristics.h"
#include <iostream>
#include <string>

class CRModel{
protected:
  Metaparameters* metaparameters;
  ntensor* eq_vals; // we allow the possibility of multiple equilibria (hence vector of vector)
  Model_parameters* model_param;
  func_equ_evol equations_of_evolution;

public:
  /* CONSTRUCTORS */
  CRModel();
  CRModel(const CRModel&);
  CRModel(Metaparameters&);
  CRModel(const foodmatrix&, Metaparameters&);
  CRModel(Model_parameters*);

  /* DESTRUCTOR */
  virtual ~CRModel();

  /* OPERATORS OVERLOAD */
  CRModel& operator=(const CRModel&);


  void create_model_parameters(Metaparameters&);
  nctype largest_eigenvalue_at_equilibrium() const;
  ncvector eigenvalues_at_equilibrium() const;
  void save(std::ostream&) const; // outputs the model to the external file
  std::ostream& display(std::ostream&) const;
  bool energy_constraint() const;
  bool constraints_fulfilled(const Metaparameters& m) const;
  bool positive_parameters() const;
  bool is_dynamically_stable() const;
  void save_simulation() const;
  void save_jacobian_at_equilibrium(std::string) const;
  void write_death_rates(std::string) const;
  Dynamical_variables perturb_equilibrium() const;
  void perturb_parameters() const;
  void perturb_parameters(const ntype &) const;
  void save_new_equilibrium(const Extinction&) const;
  bool respects_equations_of_evolution_at_equilibrium() const;

  /* virtual functions */
  virtual void attempt_to_build_model(const foodmatrix&,Metaparameters&, unsigned int);
  virtual void attempt_to_build_model_with_m(const foodmatrix&, Metaparameters&, unsigned int);
  virtual nmatrix jacobian_at_equilibrium() const;
  virtual nmatrix jacobian(const Dynamical_variables&) const; // returns the jacobian for the given dynamical variables




  Metaparameters* get_metaparameters() const;
  Model_parameters* get_model_parameters() const;
  ntensor* get_equilibrium_abundances() const;
  Parameter_set* get_parameter_set() const;
  func_equ_evol get_equations_of_evolution() const;

  double get_m0() const;
  double get_d0() const;
  nvector get_m() const;
  nvector get_d() const;

  /* returns the n-th equilibrium value of the resources */
  nvector get_resources_equilibrium(unsigned int n=0) const;
  /* returns the n-th equilibrium value of the species */
  nvector get_consumers_equilibrium(unsigned int n=0) const;

  ntype get_resilience_jacobian() const;
  ntype get_resilience_dynamical_stability(const ntype& delta=0.);


  /* gives back the Beta and Gamma matrices from the jacobian at equilibrium */
  nmatrix get_Beta_matrix(unsigned int eq_number=0) const;
  nmatrix get_Gamma_matrix(unsigned int eq_number=0) const;

  /* FLUXES PART */
  /* FLUXES FOR RESOURCES */
  /* environmental input biomass flux for resource mu */
  ntype environmental_flux_resource(unsigned int mu) const;
  ntype environmental_flux_equilibrium_resource(unsigned int mu, unsigned int equilibrium_number=0) const;
  /* diffusion flux for resource mu */
  ntype diffusion_flux_resource(unsigned int mu, const nvector& R) const ;
  ntype diffusion_flux_equilibrium_resource(unsigned int mu, unsigned int equilibrium_number=0) const;
  /* syntrophy flux for resource mu */
  ntype syntrophy_flux_resource(unsigned int mu, const nvector& R, const nvector& S) const;
  ntype syntrophy_flux_equilibrium_resource(unsigned int mu, unsigned int equilibrium_number=0) const;
  /* consumption from consumers flux for resource mu */
  ntype consumption_flux_resource(unsigned int mu, const nvector& R, const nvector& S) const;
  ntype consumption_flux_equilibrium_resource(unsigned int mu, unsigned int equilibrium_number=0) const;

  /* FLUXES FOR CONSUMERS */
  /* consumption intake flux for consumer */
  ntype consumption_intake_flux_consumer(unsigned int i, const nvector& R, const nvector& S) const;
  ntype consumption_intake_flux_equilibrium_consumer(unsigned int i, unsigned int equilibrium_number=0) const;
  /* diffusion flux for consumer */
  ntype diffusion_flux_consumer(unsigned int i, const nvector & S) const;
  ntype diffusion_flux_equilibrium_consumer(unsigned int i, unsigned int equilibrium_number=0) const;
  /* syntrophy flux for consumer */
  ntype syntrophy_flux_consumer(unsigned int i, const nvector& R, const nvector & S) const;
  ntype syntrophy_flux_equilibrium_consumer(unsigned int i, unsigned int equilibrium_number=0) const;
  /* gives back the order parameter */
  ntype order_parameter() const;

  /* tells you if the system is stable, marginally stable or unstable dynamically */
  systemstability assess_dynamical_stability() const;

  bool has_linearly_stable_eq() const;
  ntype critical_radius() const;
  bool is_in_weak_LRI() const;
  bool is_in_strong_LRI() const;

  nmatrix get_first_equilibrium() const;

  /* returns dynamical variables perturbed away from their equilibrium value*/
  nmatrix perturb_abundances(const ntype& );

  /* returns the extinction properties with the initial values of abundances */
  Extinction evolve_until_equilibrium_from_abundances(const nmatrix& , ntype threshold=1e-9, eqmode eq_mode = oneextinct, writemode write_mode=writemode()) const;
  /* returns the evolution from equilibrium (assumes we are not at equilibrium with R^* and S^*) -> finds the new eq*/
  Extinction evolve_until_equilibrium(ntype, eqmode eq_mode=convergence, writemode write_mode=writemode()) const;
  /* returns the general extinction properties for the initial values init_val */
  Extinction evolve_until_equilibrium_general(const nmatrix& init_val, ntype threshold, eqmode eq_mode, writemode write_mode=writemode())const;
};



#endif
