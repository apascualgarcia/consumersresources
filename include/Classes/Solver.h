#ifndef SOLVER_H
#define SOLVER_H

#include "Classes/Custom_types.h"
#include "Classes/Metaparameters.h"

struct Solver_Parameters{
  // these dictate the values we will take for compute_average_extinction
  Metaparameters* metaparameters;
  unsigned int Nsimul;
  eqmode equilibrium;
};

struct Delta_critical{
  ntype delta_crit;
  ntype delta_low;
  ntype delta_high;
  ntype accuracy;
};

#endif
