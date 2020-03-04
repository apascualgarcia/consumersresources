#ifndef DYNAMICAL_VARIABLES_H
#define DYNAMICAL_VARIABLES_H
#include "Classes/Custom_types.h"
class Dynamical_variables{
private :
  nvector* resources;
  nvector* consumers;
public:
  Dynamical_variables(nvector* resources_, nvector* consumers_);
  nvector* get_resources() const;
  nvector* get_consumers() const;
};
#endif
