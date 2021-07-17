#include "CRModel.h"

EcologicalNetwork::EcologicalNetwork(){
  return;
}

EcologicalNetwork::EcologicalNetwork(const unsigned int NR_, const unsigned int NS_, const double & g_conn){
  this->NR= NR_;
  this->NS = NS_;

  /* take by default A as a zero matrix first */
  this->A = nmatrix(this->NR, nvector(this->NS, 0.));
  this->G = random_binary_matrix_with_connectance(this->NS, this->NR, g_conn);
  return;
}

EcologicalNetwork::EcologicalNetwork(const Metaparameters& m){
  this->G = load_food_matrix(m);
  this->NS = this->G.size();
  this->NR = this->G[0].size();
  this->A = nmatrix(this->NR, nvector(this->NS, 0.));
  return;
}


void EcologicalNetwork::optimize(MonteCarloSolver& mcs){
  if(mcs.mcmode==constant_connectance){
    this->A = random_binary_matrix_with_connectance(this->NR, this->NS, connectance(this->G));
  }
  if(mcs.mcmode==both_modified){
    this->A = random_binary_matrix_with_connectance(this->NR, this->NS, 0.5);
    this->G = random_binary_matrix_with_connectance(this->NS, this->NR, connectance(this->G));
  }
  apply_MC_algorithm(*this, mcs);
  return;
}
