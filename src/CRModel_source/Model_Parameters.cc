#include "../../include/CRModel.h"

Model_parameters::Model_parameters(){
  params = new Parameter_set();
  return;
}
Model_parameters::~Model_parameters(){
  delete this->params;
  return;
}

Parameter_set* Model_parameters::get_parameters() const{
  return params;
}

void Model_parameters::set_sigma(const nmatrix& s) const {
  params->sigma = s;
  return;
}
void Model_parameters::set_alpha(const nmatrix& a) const {
  params->alpha = a;
  return;
}
void Model_parameters::set_gamma(const nmatrix& g) const {
  params->gamma = g;
  return;
}
void Model_parameters::set_tau(const nmatrix& t) const {
  params->tau = t;
  return;
}
void Model_parameters::set_l(const nvector& el) const {
  params->l = el;
  return;
}
void Model_parameters::set_d(const nvector& de) const {
  params->d = de;
  return;
}
void Model_parameters::set_m(const nvector& em) const {
  params->m = em;
  return;
}
void Model_parameters::set_NR(const unsigned int & enar) const {
  params->NR = enar;
  return;
}
void Model_parameters::set_NS(const unsigned int & enes) const {
  params->NS = enes;
  return;
}
