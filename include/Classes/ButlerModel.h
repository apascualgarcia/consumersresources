#ifndef BUTLER_MODEL_H
#define BUTLER_MODEL_H

#include "Consumer_Resource_Model.h"

class ButlerModel:public CRModel{
public:
  ButlerModel();
  ButlerModel(Metaparameters &);
  ButlerModel(const foodmatrix&, Metaparameters &);
  ButlerModel(const CRModel &);

  virtual void attempt_to_build_model(const foodmatrix&,Metaparameters&, unsigned int);
};


#endif
