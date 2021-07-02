#ifndef METAPARAMETERS_H
#define METAPARAMETERS_H

#include "Classes/Custom_types.h"
#include <string>

struct Metaparameters{
  ntype gamma0;
  ntype alpha0;
  ntype sigma0;
  ntype p;
  ntype R0;
  ntype S0;
  ntype l0;
  ntype m0;
  ntype epsilon;
  unsigned int NR;
  unsigned int NS;
  gammamode gamma_mode;
  taumode tau_mode;
  alphamode alpha_mode;
  std::string foodmatrixpath;
  std::string syntrophy_matrix_path;
  unsigned int verbose;
  bool energy_constraint;
  bool budget_constraint;
  unsigned int nb_attempts;
  unsigned int seed_number;
  ntype tf;
  std::string save_path;
  ntype perturb_eq;
  ntype perturb_parameters;
  eqmode equilibrium;
  ntype convergence_threshold;
  buildingmode building_mode;
  /* path to which point should be studied (eg common feasible volume) */
  std::string volume_of_interest_path;
  /* tells in which way the system should be structurally perturbed */
  unsigned int struct_pert_type;
  /* next step matrix mode for possible MC solver */
  MCmode mcmode;

  Metaparameters(int argc, char *argv[]);
  /* gives back the hard limit over which we know we won't find any feasible system*/
  ntype physical_maximum_alpha0() const;

  /* gives back a limit below which we know we will have a feasible system */
  ntype minimum_S0_guaranteed_feasability() const;


  /*  gives back the softer limit after which prob(draw feasible system) < 1,
      has an accuracy on alpha of roughly alpha_accuracy */
  ntype feasible_alpha_max(ntype alpha_accuracy = 1e-7) const;

  /*  gives back the critical dynamical syntrophy for a set of metaparameters, ie
      the largest syntrophy for which we have systems that are fully dynamically stable */
  ntype dynamical_alpha_max(ntype alpha_accuracy = 1e-7) const;
  /*  gives back the softer limit before which prob(draw feasible system) < 1,
      has an accuracy on alpha of roughly alpha_accuracy */
  ntype feasible_alpha_min(ntype alpha_accuracy = 1e-7) const;

  /* returns the maximum feasible S0 taking into account the other metaparameters */
  ntype feasible_S0_max(ntype S0_accuracy=1e-7) const;

  /* returns the maximum feasible gamma0 taking into account the other metaparameters */
  ntype feasible_gamma0_max(ntype gamma0_accuracy=1e-7)const;

  ntype quadratic_form_low_intra_resource_interaction(const nmatrix&, const nmatrix&) const;

  /* gives you back the common feasible volume : works now only for alpha0 = 0 */
  nmatrix common_feasible_volume(unsigned int Npoints) const;

  /* returns the set of (gamma0,S0) which are fully locally dynamically stable */
  nmatrix set_of_lds_points() const;
  /* returns the set of (gamma0,S0) which are fully feasible */
  nmatrix set_of_feasible_points() const;

  /* loads the points given in the path */
  nmatrix load_volume() const;

  /* attempt at building a more accurate LRI quadratic form */
  ntype accurate_quadratic_form_LRI(const nmatrix& A, const nmatrix& G) const;

  /* new attempt at building an LRI quadratic form */
  ntype newly_corrected_quadratic_form_LRI(const nmatrix& A, const nmatrix& G) const;

  /* estimates critical radius with A and G matrix as input */
  ntype critical_radius(const nmatrix& A, const nmatrix& G) const;
};


#endif
