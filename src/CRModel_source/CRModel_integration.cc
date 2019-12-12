#include "CRModel.h"
#include <iostream>
#include <fstream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <assert.h>
#include <vector>
#include <array>
#include <math.h>


nmatrix CRModel::time_evolution(const Dynamical_variables& init_val, ntype tf) const{
  // evolution contains the state of the system at all times i.e. evolution[0] contains
  // the initial value of the Dynamical_variables init_val; The goal of the algorithm is
  // to integrate until time t.
  nmatrix evolution;
  Metaparameters* p = this->metaparameters;

  // initialization of the system
  double t0=0.;
  double t = t0;
  double y[p->NR+p->NS];
  double h = 1e-6;

  for(size_t nu = 0; nu < p->NR; ++nu){
    y[nu]= (*init_val.get_resources())[nu];
  }
  for(size_t i = 0; i < p->NS; ++i){
    y[i+p->NR]=(*init_val.get_consumers())[i];
  }

  std::cout << "Initial value of system: ";
  for(size_t i = 0; i < p->NR+p->NS; ++i){
    std::cout << y[i] << " ";
  }
  std::cout << std::endl;

  // initialization of the integrator
  size_t dim = p->NR+p->NS;
  const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rkf45;
  gsl_odeiv2_step* s = gsl_odeiv2_step_alloc(T, dim);
  gsl_odeiv2_control* c = gsl_odeiv2_control_y_new(1e-6, 0.0);
  gsl_odeiv2_evolve* e = gsl_odeiv2_evolve_alloc(dim);
  gsl_odeiv2_system sys = {ode_equations_of_evolution, NULL, dim, this->model_param->get_parameters()};


  // system temporal evolution
  while(t<tf){
    int status = gsl_odeiv2_evolve_apply(e,c,s, &sys, &t, tf, &h, y);
    if (status != GSL_SUCCESS){
      std::cerr << "Error in the integration of the ODE!" << std::endl;
      break;
    }
    nvector v;
    v.push_back(t);
    for(size_t nu=0; nu < p->NR; ++nu){
      v.push_back(y[nu]);
    }
    for(size_t i=0; i < p->NS; ++i){
      v.push_back(y[i+p->NR]);
    }
    evolution.push_back(v);
  }

  gsl_odeiv2_evolve_free(e);
  gsl_odeiv2_control_free(c);
  gsl_odeiv2_step_free(s);

  return evolution;
}

// assume y contains resources then
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
void CRModel::write_time_evolution_until_equilibrium(const Dynamical_variables& init_val, ntype tmax, ntype threshold) const{
  Metaparameters* p = this->metaparameters;

  // initialization of the system
  double t0=0.;
  double t = t0;
  double y[p->NR+p->NS];
  double h = 1e-6;

  for(size_t nu = 0; nu < p->NR; ++nu){
    y[nu]= (*init_val.get_resources())[nu];
  }
  for(size_t i = 0; i < p->NS; ++i){
    y[i+p->NR]=(*init_val.get_consumers())[i];
  }

  // initialization of the integrator
  size_t dim = p->NR+p->NS;
  const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rkf45;
  gsl_odeiv2_step* s = gsl_odeiv2_step_alloc(T, dim);
  gsl_odeiv2_control* c = gsl_odeiv2_control_y_new(1e-6, 0.0);
  gsl_odeiv2_evolve* e = gsl_odeiv2_evolve_alloc(dim);
  gsl_odeiv2_system sys = {ode_equations_of_evolution, NULL, dim, this->model_param->get_parameters()};

  unsigned int counts = 0;
  double previous_y[10][p->NR+p->NS];
  double eq_coeff = 1.;

  //opening the writing file
  std::ofstream myfile;
  myfile.open(p->save_path, std::ios::app);
  bool save_success(false);
  if(not(myfile.is_open())){
    std::cerr << "Could not open " << p->save_path << " to write the temporal evolution" << std::endl;
  }else{
    save_success = true;
  }

  ncvector eigvals = this->eigenvalues_at_equilibrium();
  ntype maxeigval = real(eigvals[0]);
  for(size_t i = 1; i < eigvals.size(); ++i){
    if((real(eigvals[i]) > maxeigval) and abs(real(eigvals[i])) > EIGENSOLVER_PRECISION){
      maxeigval = real(eigvals[i]);
    }
  }

  std::cout << maxeigval << std::endl;
  myfile << "# " << maxeigval << std::endl;

  // system temporal evolution
  while( (t<tmax) and (eq_coeff > threshold)){
    int status = gsl_odeiv2_evolve_apply(e,c,s, &sys, &t, tmax, &h, y);
    if (status != GSL_SUCCESS){
      std::cerr << "Error in the integration of the ODE!" << std::endl;
      break;
    }

    // store the previous ten values to estimate convergence
    for(size_t i = 0; i < p->NR+p->NS; ++i){
      previous_y[counts%10][i] = y[i];
    }

    // stores the mean value over the last ten runs of every resource and consumer
    if(counts >= 10){
      double mean_el[p->NR+p->NS];
      for(size_t i = 0; i < p->NR+p->NS; ++i){
        mean_el[i] = 0.;
        for(size_t j = 0; j < 10; ++j){
          mean_el[i] += previous_y[j][i]/10.;
        }
      }

      // computes the "equilibrium coefficient" to estimate
      eq_coeff = 0.;
      for(size_t i = 0 ; i < p->NR+p->NS; ++i){
        eq_coeff += pow(mean_el[i]-y[i], 2.);
      }
      eq_coeff = pow(eq_coeff, 0.5);
    }
    counts += 1;
    myfile << t <<" ";
    for(size_t i = 0; i < p->NR+p->NS; ++i){
      myfile << y[i] << " ";
    }
    myfile << eq_coeff << std::endl;

  }

  myfile.close();

  gsl_odeiv2_evolve_free(e);
  gsl_odeiv2_control_free(c);
  gsl_odeiv2_step_free(s);

  return;
}

Extinction CRModel::evolve_until_equilibrium(ntype threshold, eqmode eq_mode) const{
  Metaparameters* p = this->metaparameters;
  // initialization of the system
  double t0=0.;
  double t = t0;
  double y[p->NR+p->NS];
  double h = 1e-6;
  double tmax = double(p->tf);

  for(size_t nu = 0; nu < p->NR; ++nu){
    y[nu]= (*eq_vals)[0][0][nu];
  }
  for(size_t i = 0; i < p->NS; ++i){
    y[i+p->NR]=(*eq_vals)[0][1][i];
  }

  // initialization of the integrator
  size_t dim = p->NR+p->NS;
  const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rkf45;
  gsl_odeiv2_step* s = gsl_odeiv2_step_alloc(T, dim);
  gsl_odeiv2_control* c = gsl_odeiv2_control_y_new(1e-6, 0.0);
  gsl_odeiv2_evolve* e = gsl_odeiv2_evolve_alloc(dim);
  gsl_odeiv2_system sys = {ode_equations_of_evolution, NULL, dim, this->model_param->get_parameters()};

  unsigned int counts = 0;
  double previous_y[10][p->NR+p->NS];
  double eq_coeff = 1.;
  bool exit_condition = true;
  bool one_extinct = false;

  exit_condition = (t>=tmax) or (eq_coeff <= threshold);

  // system temporal evolution
  while(not(exit_condition)){
    int status = gsl_odeiv2_evolve_apply(e,c,s, &sys, &t, tmax, &h, y);
    if (status != GSL_SUCCESS){
      std::cerr << "Error in the integration of the ODE!" << std::endl;
      break;
    }

    if(p->verbose > 3){
      std::cout << "Population of the system at time t = "<< t << " :" << std::endl;
      std::cout << "  Resources :";
      for(size_t nu = 0; nu < p->NR; ++nu){
        std::cout << " " << y[nu] ;
      }
      std::cout << std::endl;
      std::cout << "  Consumers :";
      for(size_t i = p->NR; i < p->NR+p->NS; ++i){
        std::cout << " " << y[i];
      }
      std::cout << std::endl;
    }

    // store the previous ten values to estimate convergence
    for(size_t i = 0; i < p->NR+p->NS; ++i){
      previous_y[counts%10][i] = y[i];
    }

    // stores the mean value over the last ten runs of every resource and consumer
    if(counts >= 10){
      double mean_el[p->NR+p->NS];
      for(size_t i = 0; i < p->NR+p->NS; ++i){
        mean_el[i] = 0.;
        for(size_t j = 0; j < 10; ++j){
          mean_el[i] += previous_y[j][i]/10.;
        }
      }

      // computes the "equilibrium coefficient" to estimate
      eq_coeff = 0.;
      for(size_t i = 0; i < p->NR+p->NS; ++i){
        if(y[i] > 0){
          eq_coeff += pow(mean_el[i]/y[i]-1., 2.);
        }
      }
      eq_coeff = pow(eq_coeff, 0.5)/(p->NR+p->NS);
    }
    counts += 1;

    /* if a resource/consumer is too small, we effectively set it to zero */
    for(size_t i=0; i < p->NR+p->NS and not(one_extinct); ++i){
      if(y[i] < threshold){
        y[i]=0.;
        if(i >= p->NR){
          one_extinct=true;
          if(p->verbose > 2){
            std::cout << " Consumer " << i-p->NR << " went extinct at time t=" << t << std::endl;
          }
        }else{
          if(p->verbose > 2){
            std::cout << " Species " << i << " went extinct at time t=" << t << std::endl;
          }
        }
        }
    }

    exit_condition = (t>=tmax) or (eq_coeff <= threshold);
    if(eq_mode == oneextinct){
      exit_condition = exit_condition or one_extinct;
    }
  }

  if(this->metaparameters->verbose>1){
    if(eq_mode == oneextinct){
      std::cout << " Time to get to equilibrium : " << t ;
      std::cout << " (";
      if(one_extinct){
        std::cout << "observed at least one extinction";
      }else if(t >= tmax){
        std::cout << "reached max time integration";
      }else if(eq_coeff <= threshold){
        std::cout << "convergence observed";
      }
      std::cout << ")" << std::endl;
    }else if(eq_mode == convergence){
      std::cout << " Time to get convergence on all the resources/consumers time lines : " << t << std::endl;
    }
  }

  /*  At the end of the integration, we compute how many species went extinct
      This should be one if eq_mode is oneextinct (except if two species or more
      are deemed extinct at the same time) */
  unsigned int extinct(0);
  for(size_t i = p->NR; i < p->NR+p->NS; ++i){
    if(y[i] < threshold){
      extinct += 1;
    }
  }

  nvector new_R, new_S;
  for(size_t nu = 0; nu < p->NR; ++nu){
    new_R.push_back(y[nu]);
  }
  for(size_t i = p->NR; i < p->NR+p->NS; ++i){
    new_S.push_back(y[i]);
  }

  gsl_odeiv2_evolve_free(e);
  gsl_odeiv2_control_free(c);
  gsl_odeiv2_step_free(s);

  Extinction to_return = {t,ntype(extinct), new_R, new_S};
  return to_return;
}

Extinction CRModel::evolve_until_equilibrium_from_abundances(const nmatrix& init_val, ntype threshold, eqmode eq_mode) const{
  Metaparameters* p = this->metaparameters;
  /* initialization of the system */
  double t0=0.;
  double t = t0;
  double y[p->NR+p->NS];
  double h = 1e-6;
  double tmax = double(p->tf);

  nmatrix old_eq = this->get_first_equilibrium();

  for(size_t nu = 0; nu < p->NR; ++nu){
    y[nu]= init_val[0][nu];
  }
  for(size_t i = 0; i < p->NS; ++i){
    y[i+p->NR]=init_val[1][i];
  }

  /* initialization of the integrator */
  size_t dim = p->NR+p->NS;
  const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rkf45;
  gsl_odeiv2_step* s = gsl_odeiv2_step_alloc(T, dim);
  gsl_odeiv2_control* c = gsl_odeiv2_control_y_new(1e-6, 0.0);
  gsl_odeiv2_evolve* e = gsl_odeiv2_evolve_alloc(dim);
  gsl_odeiv2_system sys = {ode_equations_of_evolution, NULL, dim, this->model_param->get_parameters()};

  unsigned int counts = 0;
  double previous_y[10][p->NR+p->NS];
  double eq_coeff = 1.;
  bool exit_condition = true;
  bool one_extinct = false;

  exit_condition = (t>=tmax) or (eq_coeff <= threshold);

  // system temporal evolution
  while(not(exit_condition)){
    int status = gsl_odeiv2_evolve_apply(e,c,s, &sys, &t, tmax, &h, y);
    if (status != GSL_SUCCESS){
      std::cerr << "Error in the integration of the ODE!" << std::endl;
      break;
    }

    if(p->verbose > 3){
      std::cout << "Population of the system at time t = "<< t << " :" << std::endl;
      std::cout << "  Resources :";
      for(size_t nu = 0; nu < p->NR; ++nu){
        std::cout << " " << y[nu] ;
      }
      std::cout << std::endl;
      std::cout << "  Consumers :";
      for(size_t i = p->NR; i < p->NR+p->NS; ++i){
        std::cout << " " << y[i];
      }
      std::cout << std::endl;
    }

    // store the previous ten values to estimate convergence
    for(size_t i = 0; i < p->NR+p->NS; ++i){
      previous_y[counts%10][i] = y[i];
    }

    // stores the mean value over the last ten runs of every resource and consumer
    if(counts >= 10){
      double mean_el[p->NR+p->NS];
      for(size_t i = 0; i < p->NR+p->NS; ++i){
        mean_el[i] = 0.;
        for(size_t j = 0; j < 10; ++j){
          mean_el[i] += previous_y[j][i]/10.;
        }
      }

      // computes the "equilibrium coefficient" to estimate
      eq_coeff = 0.;
      for(size_t i = 0; i < p->NR+p->NS; ++i){
        if(y[i] > 0){
          eq_coeff += pow(mean_el[i]/y[i]-1., 2.);
        }
      }
      eq_coeff = pow(eq_coeff, 0.5)/(p->NR+p->NS);
    }
    counts += 1;

    /* if a resource/consumer is too small, we effectively set it to zero */
    for(size_t i=0; i < p->NR+p->NS and not(one_extinct); ++i){
      if(y[i] < threshold){
        y[i]=0.;
        if(i >= p->NR){
          one_extinct=true;
          if(p->verbose > 2){
            std::cout << " Consumer " << i-p->NR << " went extinct at time t=" << t << std::endl;
          }
        }else{
          if(p->verbose > 2){
            std::cout << " Species " << i << " went extinct at time t=" << t << std::endl;
          }
        }
        }
    }

    exit_condition = (t>=tmax) or (eq_coeff <= threshold);
    if(eq_mode == oneextinct){
      exit_condition = exit_condition or one_extinct;
    }
  }

  if(this->metaparameters->verbose>1){
    if(eq_mode == oneextinct){
      std::cout << " Time to get to equilibrium : " << t ;
      std::cout << " (";
      if(one_extinct){
        std::cout << "observed at least one extinction";
      }else if(t >= tmax){
        std::cout << "reached max time integration";
      }else if(eq_coeff <= threshold){
        std::cout << "convergence observed";
      }
      std::cout << ")" << std::endl;
    }else if(eq_mode == convergence){
      std::cout << " Time to get convergence on all the resources/consumers time lines : " << t << std::endl;
    }
  }

  /*  At the end of the integration, we compute how many species went extinct
      This should be one if eq_mode is oneextinct (except if two species or more
      are deemed extinct at the same time) */
  unsigned int extinct(0);
  for(size_t i = p->NR; i < p->NR+p->NS; ++i){
    if(y[i] < threshold){
      extinct += 1;
    }
  }

  nvector new_R, new_S;
  for(size_t nu = 0; nu < p->NR; ++nu){
    new_R.push_back(y[nu]);
  }
  for(size_t i = p->NR; i < p->NR+p->NS; ++i){
    new_S.push_back(y[i]);
  }

  gsl_odeiv2_evolve_free(e);
  gsl_odeiv2_control_free(c);
  gsl_odeiv2_step_free(s);

  Extinction to_return = {t,ntype(extinct), new_R, new_S, old_eq[0], old_eq[1]};
  return to_return;

}
