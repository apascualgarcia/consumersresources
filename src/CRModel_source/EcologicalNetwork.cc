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

nmatrix EcologicalNetwork::get_effective_competition_matrix(const Metaparameters& m) const{
  nmatrix C = nmatrix(m.NS, nvector(m.NS, 0.)), G=this->G, A=this->A;
  nvector D = nvector(m.NR, 0.);

  /* initialize D */
  for(size_t nu=0; nu < m.NR; ++nu){
    for(size_t k=0; k < m.NS; ++k){
      D[nu]+=A[nu][k];
    }
    D[nu]*=(1.*m.alpha0*m.S0);
    D[nu]+=m.l0;
    D[nu]/=(1.*m.R0);
  }

  for(size_t i=0; i < m.NS; ++i){
    for(size_t j=0; j < m.NS; ++j){
      for(size_t nu=0; nu < m.NR; ++nu){
        ntype G_sum=0.;
        for(size_t k=0; k < m.NS; ++k){
          G_sum+=G[k][nu];
        }
        C[i][j]+=m.sigma0*m.gamma0/(D[nu]*D[nu])*G[i][nu]*(m.gamma0*m.l0*G[j][nu]-m.alpha0*A[nu][j]*(D[nu]+m.gamma0*m.S0*G_sum));
      }
    }
  }

  return C;
}

ntype EcologicalNetwork::effective_competition(const Metaparameters& m) const{
  nmatrix C = this->get_effective_competition_matrix(m);
  ntype lambda1 = real(largest_eigenvalue(C));

  return 1./(m.NS-1.)*((1.*m.NS*lambda1)/trace(C)-1.);
}
