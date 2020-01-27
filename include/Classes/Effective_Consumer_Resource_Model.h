#ifndef EFFECTIVE_CONSUMER_RESOURCE_MODEL_H
#define EFFECTIVE_CONSUMER_RESOURCE_MODEL_H

#include "Consumer_Resource_Model.h"

class EffectiveCRModel:public CRModel{
public:
  EffectiveCRModel(Metaparameters&);
  EffectiveCRModel(const foodmatrix&, Metaparameters&);
  EffectiveCRModel(const CRModel &);
};

#endif
