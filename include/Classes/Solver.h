#ifndef SOLVER_H
#define SOLVER_H

#include "Classes/Custom_types.h"
#include "Classes/Metaparameters.h"
#include <gsl/gsl_math.h>

struct Solver_Parameters{
  // these dictate the values we will take for compute_average_extinction
  Metaparameters* metaparameters;
  unsigned int Nsimul;
  eqmode equilibrium;
  stabilitymode stab_mode; // either structural or dynamical
  double target; //basically allows us to say we want to solve f(x)=shift instead of f(x)=0
};

struct Delta_critical{
  ntype delta_crit;
  ntype delta_low;
  ntype delta_high;
  ntype accuracy;
};

struct delta_solver{
  fitmode fit_mode;
  eqmode eq_mode;
  stabilitymode stab_mode;
};


#endif
