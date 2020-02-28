#include "../../include/CRModel.h"

ButlerModel::ButlerModel(){
  this->metaparameters = NULL;
  this->eq_vals = NULL;
  this->model_param = NULL;
  this->equations_of_evolution = ode_equations_of_evolution;
}
ButlerModel::ButlerModel(Metaparameters& meta){
  this->equations_of_evolution=ode_equations_of_evolution;
  unsigned int attempts(0);
  this->create_model_parameters(meta);
  foodmatrix F=load_food_matrix(meta);
  do{
    attempts+=1;
    this->attempt_to_build_model(F, meta, attempts);
  }while(not(this->constraints_fulfilled(meta)));

  if(meta.verbose > 1){
    std::cout << "\t Feasible system built in "<<attempts<<" iteration(s). ";
    if(attempts > meta.nb_attempts){
      std::cout << "\t The metaparameters had to be changed.";
    }else{
      std::cout << "It was possible to use the initial metaparameters.";
    }
    std::cout << std::endl;
  }

  return;
}

ButlerModel::ButlerModel(const foodmatrix& F, Metaparameters& meta){
  this->equations_of_evolution=ode_equations_of_evolution;
  unsigned int attempts(0);
  this->create_model_parameters(meta);
  do{
    attempts+=1;
    this->attempt_to_build_model(F, meta, attempts);
  }while(not(this->constraints_fulfilled(meta)));

  if(meta.verbose > 1){
    std::cout << "\t Feasible system built in "<<attempts<<" iteration(s). ";
    if(attempts > meta.nb_attempts){
      std::cout << "\t The metaparameters had to be changed.";
    }else{
      std::cout << "It was possible to use the initial metaparameters.";
    }
    std::cout << std::endl;
  }
  return;
}
ButlerModel::ButlerModel(const CRModel& model){
  this->metaparameters= &(*(model.get_metaparameters()));
  this->eq_vals = new ntensor(*(model.get_equilibrium_abundances()));
  this->model_param = new Model_parameters(*(model.get_model_parameters()));
  this->equations_of_evolution = func_equ_evol(model.get_equations_of_evolution());
}


void ButlerModel::attempt_to_build_model(const foodmatrix& F, Metaparameters& meta, unsigned int attempts){
  Parameter_set* p = this->model_param->get_parameters();

  /* first find the values for the equilibria */
  nvector Req = build_resources(meta);
  nvector Seq = build_consumers(meta);

  nmatrix equilibria;
  equilibria.push_back(Req);
  equilibria.push_back(Seq);

  ntensor equilibria_vals;
  equilibria_vals.push_back(equilibria);

  *(this->eq_vals) = equilibria_vals;

  p->NR = meta.NR;
  p->NS = meta.NS;

  /* first sigma, Req, Seq are drawn randomly */
  p->sigma = build_sigma_Butler(meta);

  /* then we build gamma according to the food matrix */
  p->gamma = build_gamma(F,meta);

  /* alpha and tau are set to zero */
  p->alpha = build_alpha(p, meta, Req, attempts);
  p->tau = build_tau(p, meta, attempts);

  /* d is then set */
  nvector d;
  for (size_t i=0; i < meta.NS; ++i){
    ntype result = 0.;
    for (size_t mu =0 ; mu < meta.NR; ++mu){
      result+=(p->sigma)[i][mu]*(p->gamma)[i][mu]*Req[mu]-(p->tau)[mu][i];
    }
    d.push_back(result);
  }
  p->d = d;

  /* with Butler's implementation, since m=alpha=0, we do not have any freedom over l anymore */
  nvector l, m(meta.NR, 0.);
  for(size_t mu=0; mu < meta.NR; ++mu){
    ntype result = 0.;
    for(size_t j=0; j < meta.NS; ++j){
      result+=(p->gamma[j][mu])*(Req[mu])*(Seq[j]);
    }

    l.push_back(result);
  }

  p->l = l;
  p->m = m;

  return;
}
