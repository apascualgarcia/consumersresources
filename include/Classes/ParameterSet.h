#ifndef PARAMETER_SET_H
#define PARAMETER_SET_H
#include "Classes/Custom_types.h"

struct Parameter_set{
  nmatrix alpha;
  nmatrix gamma;
  nmatrix tau;
  nmatrix sigma;
  nvector l;
  nvector m;
  nvector d;
  unsigned int NR;
  unsigned int NS;
};

#endif
