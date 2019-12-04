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
  ntype epsilon;
  unsigned int NR;
  unsigned int NS;
  gammamode gamma_mode;
  taumode tau_mode;
  alphamode alpha_mode;
  std::string foodmatrixpath;
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

  Metaparameters(int argc, char *argv[]);
};

#endif
