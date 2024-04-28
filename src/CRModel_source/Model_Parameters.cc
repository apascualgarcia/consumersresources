#include "../../include/CRModel.h"

Model_parameters::Model_parameters(){
  params = Parameter_set();
  return;
}

Model_parameters::Model_parameters(const Model_parameters& other_mod){
  this->params = other_mod.get_parameter_set();
  return;
}

Model_parameters::~Model_parameters(){
  //delete this->params;
  return;
}

Model_parameters::Model_parameters(const Metaparameters& meta, unsigned int attempts){

  this->params.NR = meta.NR;
  this->params.NS = meta.NS;

  /* first sigma, Req, Seq are drawn randomly */
  this->params.sigma = build_sigma_Butler(meta);

  /* then we build gamma according to the food matrix */
  this->params.gamma = build_gamma(load_food_matrix(meta),meta);

  /* then we build alpha according to the other parameters */
  /* first find the values for the equilibria */
  nvector Req = build_resources(meta);
  nvector Seq = build_consumers(meta);

  this->params.alpha = build_alpha(&(this->params), meta, Req, attempts);
  this->params.tau = build_tau(&(this->params), meta, attempts);

  /* d is then set */
  nvector d;
  for (size_t i=0; i < meta.NS; ++i){
    ntype result = 0.;
    for (size_t mu =0 ; mu < meta.NR; ++mu){
      result+=(this->params.sigma)[i][mu]*(this->params.gamma)[i][mu]*Req[mu]-(this->params.tau)[mu][i];
    }
    d.push_back(result);
  }

  this->params.d = d;

  /* still have to set l and m */
  nvector l, m(meta.NR);
  l = build_l(meta);
  for(size_t nu=0; nu < this->params.NR; ++nu){
    ntype C = 0.;
    for(size_t j = 0; j < this->params.NS; ++j){
      C+=(this->params.alpha[nu][j]*Seq[j]-this->params.gamma[j][nu]*Req[nu]*Seq[j]);
    }
    m[nu]=(ntype(l[nu]+C)/Req[nu]);
  }
  this->params.l = l;
  this->params.m = m;
  return;
}

void Model_parameters::display(std::ostream& os) const{
  os << params;
  return;
}

Parameter_set Model_parameters::get_parameter_set() const{
  return params;
}

Parameter_set* Model_parameters::get_parameters(){
  return &(this->params);
}

void Model_parameters::set_sigma(const nmatrix& s)  {
  params.sigma = s;
  return;
}
void Model_parameters::set_alpha(const nmatrix& a)  {
  params.alpha = a;
  return;
}
void Model_parameters::set_gamma(const nmatrix& g)  {
  params.gamma = g;
  return;
}
void Model_parameters::set_tau(const nmatrix& t)  {
  params.tau = t;
  return;
}
void Model_parameters::set_l(const nvector& el)  {
  params.l = el;
  return;
}
void Model_parameters::set_d(const nvector& de)  {
  params.d = de;
  return;
}
void Model_parameters::set_m(const nvector& em)  {
  params.m = em;
  return;
}
void Model_parameters::set_NR(const unsigned int & enar)  {
  params.NR = enar;
  return;
}
void Model_parameters::set_NS(const unsigned int & enes)  {
  params.NS = enes;
  return;
}


void Model_parameters::optimize(MonteCarloSolver & mcs, void* extra_params) {
  mcs.optimization_procedure(*this, extra_params);
  return;
}