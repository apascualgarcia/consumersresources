#ifndef EXTINCTION_H
#define EXTINCTION_H
#include "Classes/Custom_types.h"

struct Extinction{
  ntype t_eq;
  ntype extinct;
  nvector new_Req;
  nvector new_Seq;
};

struct Extinction_statistics{
  statistics t_eq;
  statistics extinct;
  nvector_statistics new_Req;
  nvector_statistics new_Seq;
};


#endif
