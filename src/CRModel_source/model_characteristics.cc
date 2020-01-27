#include "CRModel.h"

int ode_equations_of_evolution(double t, const double y[], double f[], void* params){
  Parameter_set* p = &(*(Parameter_set*) params);
  nvector R;
  nvector S;

  /*  we fill S and R and set every member of S and R to zero
      if it's smaller than the precision of the integrator */
  for(size_t nu=0; nu < p->NR ; ++nu){
    R.push_back(y[nu]);
    if(abs(R[nu])<INTEGRATOR_ABS_PRECISION){
      R[nu]=0.;
    }
  }

  for(size_t i=0; i < p->NS; ++i){
    S.push_back(y[i+p->NR]);
    if(abs(S[i])<INTEGRATOR_ABS_PRECISION){
      S[i]=0.;
    }
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

  for(size_t mu=0; mu < p->NR;++mu){
    R.push_back(y[mu]);
    f[mu]=0.;
  }
  for(size_t i=0; i < p->NS;++i){
    S.push_back(y[i+p->NR]);
    if(abs(S[i])<INTEGRATOR_ABS_PRECISION){
      S[i]=0.;
    }
  }

  for(size_t i=0; i < p->NS; ++i){
    double result=0.;
    /* first term */
    for(size_t nu=0; nu < p->NR; ++nu){
      double denom = p->m[nu];
      for(size_t k=0; k < p->NS; ++k){
        denom+=p->gamma[k][nu]*S[k];
      }
      result+=((p->sigma[i][nu]*p->gamma[i][nu]*p->l[nu])/denom-p->alpha[nu][i]);
    }
    result-=p->d[i];
    for(size_t j=0; j < p->NS; ++j){
      double Mij =0.;
      for(size_t nu=0; nu < p->NR;++nu){
        double denom = p->m[nu];
        for(size_t k=0; k < p->NS;++k){
          denom+=p->gamma[k][nu]*S[k];
        }
        Mij+=(p->sigma[i][nu]*p->gamma[i][nu]*p->alpha[nu][j])/denom;
      }
      result+=Mij*S[j];
    }
    result*=S[i];
    f[p->NR+i]=result;
  }

  /* we implement here the effective equations of evolution, i.e. we assume dR/dt=0 and the resource evolution is not modelled */
  /* since the resources are constant we can just ignore their evolution */



  return GSL_SUCCESS;
}
