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
