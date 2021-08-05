#include "CRModel.h"

EcologicalNetwork::EcologicalNetwork(){
  return;
}

EcologicalNetwork::EcologicalNetwork(const unsigned int NR_, const unsigned int NS_, const double & g_conn){
  this->NR= NR_;
  this->NS = NS_;

  /* take by default A as a zero matrix first */
  this->A = nmatrix(this->NR, nvector(this->NS, 0.));
  this->G = random_full_rank_binary_matrix_with_connectance(this->NS, this->NR, g_conn);
  return;
}

EcologicalNetwork::EcologicalNetwork(const Metaparameters& m){
  this->G = load_food_matrix(m);
  this->NS = this->G.size();
  this->NR = this->G[0].size();
  this->A = nmatrix(this->NR, nvector(this->NS, 0.));
  return;
}

EcologicalNetwork::EcologicalNetwork(const nmatrix& A_, const nmatrix& G_){
  this->A = A_;
  this->G = G_;

  this->NS = this->G.size();
  this->NR = this->A.size();
  return;
}


void EcologicalNetwork::optimize(MonteCarloSolver& mcs){
  if(mcs.mcmode==constant_connectance){
    this->A = random_binary_matrix_with_connectance(this->NR, this->NS, connectance(this->G));
  }
  if(mcs.mcmode==both_modified){
    this->A = random_binary_matrix_with_connectance(this->NR, this->NS, 0.5);
    this->G = random_full_rank_binary_matrix_with_connectance(this->NS, this->NR, connectance(this->G));
  }
  apply_MC_algorithm(*this, mcs);
  return;
}

ntype EcologicalNetwork::effective_competition(const Metaparameters& m) const{
  nmatrix G = this->G, A = this->A;
  ntype sigma0 = m.sigma0, R0 = m.R0, S0=m.S0, l0 = m.l0, alpha0=m.alpha0, gamma0=m.gamma0;
  nmatrix BetaGamma = sigma0*S0*gamma0*G*(-gamma0*R0*G+alpha0*A);
  nmatrix C = -R0/(S0*l0)*BetaGamma;

  return mean(C);
}
