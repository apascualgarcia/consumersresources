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
  nmatrix BetaGamma = m.sigma0*m.gamma0*m.S0*(m.alpha0*this->G*this->A-m.R0*this->G*transpose(this->G));
  ntype eff_comp=0.;
  for(size_t i =0; i < m.NS; ++i){
    for(size_t j=0; j < m.NS; ++j){
      eff_comp += BetaGamma[i][j];
    }
  }

  eff_comp*=2./(m.NS*(1.-m.NS));

  return eff_comp;
}
