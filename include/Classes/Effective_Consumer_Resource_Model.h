#ifndef EFFECTIVE_CONSUMER_RESOURCE_MODEL_H
#define EFFECTIVE_CONSUMER_RESOURCE_MODEL_H

#include "Consumer_Resource_Model.h"

class EffectiveCRModel:public CRModel{
public:
  EffectiveCRModel();
  EffectiveCRModel(Metaparameters&);
  EffectiveCRModel(const foodmatrix&, Metaparameters&);
  EffectiveCRModel(const CRModel &);

  /* implements the effective jacobian at equilibrium*/
  virtual nmatrix jacobian_at_equilibrium() const;
  virtual nmatrix jacobian(const Dynamical_variables&)const;
};

#endif
