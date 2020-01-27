#include "../../include/CRModel.h"

EffectiveCRModel::EffectiveCRModel(Metaparameters& meta):CRModel(meta){
  this->equations_of_evolution = &effective_ode_equations_of_evolution;
}
EffectiveCRModel::EffectiveCRModel(const foodmatrix& F, Metaparameters& meta):CRModel(F, meta){
  this->equations_of_evolution = &effective_ode_equations_of_evolution;
}
EffectiveCRModel::EffectiveCRModel(const CRModel & model){
  this->metaparameters = model.get_metaparameters();
  this->eq_vals=model.get_equilibrium_abundances();
  this->model_param=model.get_model_parameters();
  this->equations_of_evolution=&effective_ode_equations_of_evolution;
}
