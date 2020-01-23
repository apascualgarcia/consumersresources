#include "CRModel.h"

int ode_equations_of_evolution(double t, const double y[], double f[], void* params){
  Parameter_set* p = &(*(Parameter_set*) params);
  nvector R;
  nvector S;
  for(size_t nu=0; nu < p->NR ; ++nu){
    R.push_back(y[nu]);
  }

  for(size_t i=0; i < p->NS; ++i){
    S.push_back(y[i+p->NR]);
  }

  for (size_t nu=0; nu < p->NR; ++nu){
    ntype result(0.);
    result+=p->l[nu];
    result-=p->m[nu]*R[nu];
    for (size_t j=0; j < p->NS; ++j){
      result-=p->gamma[j][nu]*R[nu]*S[j];
      result+=p->alpha[nu][j]*S[j];
    }
    f[nu] = result;
  }

  for (size_t i=0; i<p->NS; ++i){
    ntype result=0.;
    for (size_t mu=0; mu < p->NR; ++mu){
      result+=p->sigma[i][mu]*p->gamma[i][mu]*S[i]*R[mu];
      result-=p->tau[mu][i]*S[i];
    }
    result-=p->d[i]*S[i];
    f[i+p->NR] = result;
  }

  return GSL_SUCCESS;
}


int effective_ode_equations_of_evolution(double t, const double y[], double f[], void* params){
  Parameter_set* p = &(*(Parameter_set*) params);
  nvector R;
  nvector S;

  for(size_t nu=0; nu < p->NR ; ++nu){
    R.push_back(y[nu]);
  }

  for(size_t i=0; i < p->NS; ++i){
    S.push_back(y[i+p->NR]);
  }

  /* we implement here the effective equations of evolution, i.e. we assume dR/dt=0 */
  for(size_t nu=0; nu < p->NR;++nu){
    f[nu]=0.;
  }

  for (size_t i=0; i<p->NS; ++i){
    ntype result=0.;
    for (size_t mu=0; mu < p->NR; ++mu){
      result+=p->sigma[i][mu]*p->gamma[i][mu]*S[i]*R[mu];
      result-=p->tau[mu][i]*S[i];
    }
    result-=p->d[i]*S[i];
    f[i+p->NR] = result;
  }
  return GSL_SUCCESS;
}
