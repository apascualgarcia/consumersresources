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
#include <algorithm>


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
    double mean_el[p->NR+p->NS];
    unsigned int values_to_compare;
    if(counts>=10){
      values_to_compare = 10;
    }else{
      values_to_compare=counts;
    }
    for(size_t i = 0; i < p->NR+p->NS; ++i){
      mean_el[i] = 0.;
      for(size_t j = 0; j < values_to_compare; ++j){
        mean_el[i] += previous_y[j][i]/values_to_compare;
      }
    }
    // computes the "equilibrium coefficient" to estimate
    eq_coeff = 0.;
    for(size_t i = 0 ; i < p->NR+p->NS; ++i){
      eq_coeff += pow(mean_el[i]-y[i], 2.);
    }
    eq_coeff = pow(eq_coeff, 0.5);
    std::cout << "/!\ Correct write_time_evolution_until_equilibrium" << std::endl;
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
Extinction CRModel::evolve_until_equilibrium_general(const nmatrix& init_val, ntype threshold, eqmode eq_mode) const{
  Metaparameters* p = this->metaparameters;
  /* initialization of the system */
  double t0=0.;
  double t = t0;
  double y[p->NR+p->NS];
  double step_size = 1e-6;
  double tmax = double(p->tf);

  nvector old_R = init_val[0];
  nvector old_S = init_val[1];

  for(size_t nu = 0; nu < p->NR; ++nu){
    y[nu]= old_R[nu];
  }
  for(size_t i = 0; i < p->NS; ++i){
    y[i+p->NR]=old_S[i];
  }

  /* initialization of the integrator */
  size_t dim = p->NR+p->NS;
  /* we choose an adaptative step for a rk order 4 scheme */
  const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rkf45;
  gsl_odeiv2_step* s = gsl_odeiv2_step_alloc(T, dim);
  gsl_odeiv2_control* c = gsl_odeiv2_control_y_new(INTEGRATOR_ABS_PRECISION, INTEGRATOR_REL_PRECISION);
  gsl_odeiv2_evolve* e = gsl_odeiv2_evolve_alloc(dim);
  gsl_odeiv2_system sys = {ode_equations_of_evolution, NULL, dim, this->model_param->get_parameters()};

  unsigned int counts = 0;
  double previous_y[10][p->NR+p->NS];
  double diff[10][p->NR+p->NS];
  double eq_coeff = 1.;
  bool exit_condition = false;
  bool one_extinct = false;
  /* if it turns out that the system was ALREADY at equilibrium, we will return 0 as the final time */
  bool started_at_equilibrium = false;
  std::vector<unsigned int> extinct_variables;
  double t_eq = 0.;

  exit_condition = (t>=tmax) or (eq_coeff <= threshold) or started_at_equilibrium;

  /* system temporal evolution */
  while(not(exit_condition)){

    /* store the previous ten values to estimate convergence */
    for(size_t i = 0; i < p->NR+p->NS; ++i){
      previous_y[counts%10][i] = y[i];
    }

    /* evolve the system to the next state */
    int status = gsl_odeiv2_evolve_apply(e,c,s, &sys, &t, tmax, &step_size, y);
    if (status != GSL_SUCCESS){
      std::cerr << "Error in the integration of the ODE!" << std::endl;
      abort();
      break;
    }

    /* we then compute the differences between this step and the previous one to check if we started already at equilibrium */
    if(counts <=10){
      for(size_t i=0; i < p->NR+p->NS; ++i){
        diff[counts%10][i] = y[i]-previous_y[counts%10][i];
      }
    }
    /* after ten counts we check whether the system was de facto at equilibrium */
    if(counts==10){
      started_at_equilibrium = true;
      for(size_t j=0; j < 10 and started_at_equilibrium; ++j){
        for(size_t k=0; (k < p->NR+p->NS) and started_at_equilibrium; ++k){
          double difference = abs(diff[j][k]);
          if(difference> INTEGRATOR_ABS_PRECISION){
            started_at_equilibrium = false;
          }
        }
      }
      if(started_at_equilibrium and p->verbose > 3){
        std::cout << "\t \t \t It turns out the system was already at equilibrium!" << std::endl;
      }
    }

    /* displays the values if asked */
    if(p->verbose > 3){
      std::cout << "\t \t \t Population of the system at time t = "<< t << " :" << std::endl;
      std::cout << "\t \t \t Resources :";
      for(size_t nu = 0; nu < p->NR; ++nu){
        std::cout << " " << y[nu] ;
      }
      std::cout << std::endl;
      std::cout << "\t \t \t Consumers :";
      for(size_t i = p->NR; i < p->NR+p->NS; ++i){
        std::cout << " " << y[i];
      }
      std::cout << std::endl;
    }

    /* if counts >= 10, stores the mean value over the last ten runs of every resource and consumer */
    /* this is done in order to compute the convergence coefficient */
    if(counts>=10){
      double mean_el[p->NR+p->NS];
      unsigned int values_to_compare=10;
      for(size_t i = 0; i < p->NR+p->NS; ++i){
        mean_el[i] = 0.;
        for(size_t j = 0; j < values_to_compare; ++j){
          mean_el[i] += previous_y[j][i]/values_to_compare;
        }
      }

      /* computes the "convergence coefficient" to estimate */
      eq_coeff = 0.;
      for(size_t i = 0; i < p->NR+p->NS; ++i){
        if(y[i] > 0){
          eq_coeff += pow(mean_el[i]/y[i]-1., 2.);
        }
      }
      eq_coeff = pow(eq_coeff, 0.5)/(p->NR+p->NS);
    }

    /* if a resource/consumer is too small, we effectively set it to zero */
    bool local_exit = false;
    for(size_t i=0; i < p->NR+p->NS and not(local_exit); ++i){
      if(y[i] < threshold){
        y[i]=0.;
        bool already_extinct;
        /* check if species is already extinct */
        if(std::find(extinct_variables.begin(), extinct_variables.end(), i) != extinct_variables.end()){
          /* in that case the species is already extinct */
          already_extinct = true;
        }else{
          extinct_variables.push_back(i);
          already_extinct = false;
        }
        if(i >= p->NR){
          one_extinct=true;
          if(eq_mode==oneextinct){
            local_exit=true;
          }
          if(p->verbose > 2){
            if(not(already_extinct)){
              std::cout << "\t \t Consumer " << i-p->NR << " went extinct at time t=" << t << std::endl;
            }
            if(eq_mode==oneextinct){
              std::cout << "\t \t Stopped calculating who goes extinct since we found one consumer extinct " <<std::endl;
            }
          }
        }else{
          if(p->verbose > 2){
            if(not(already_extinct)){
              std::cout << "\t \t Species " << i << " went extinct at time t=" << t << std::endl;
            }
            if(eq_mode==oneextinct){
              std::cout << "\t \t Stopped calculating who goes extinct since we found one resource extinct " <<std::endl;
            }
          }
        }
        }
    }

    /* let's not forget to add a count since we did a whole step */
    counts += 1;

    /* we check whether we should continue integrating or not */
    exit_condition = (t>=tmax) or (eq_coeff <= threshold) or started_at_equilibrium;
    if(eq_mode == oneextinct){
      exit_condition = exit_condition or one_extinct;
    }
  }

  if(not(started_at_equilibrium)){
    t_eq = t;
  }else{
    t_eq = 0.;
  }


  if(this->metaparameters->verbose>1){
    if(eq_mode == oneextinct){
      std::cout << "\t Time to get to equilibrium : " << t ;
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
      std::cout << "\t Time to get convergence on all the resources/consumers time lines : " << t << std::endl;
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

  Extinction to_return = {t_eq,ntype(extinct), new_R, new_S, old_R, old_S};
  return to_return;
}

Extinction CRModel::evolve_until_equilibrium(ntype threshold, eqmode eq_mode) const{
  return evolve_until_equilibrium_general((*eq_vals)[0],threshold, eq_mode);
}
Extinction CRModel::evolve_until_equilibrium_from_abundances(const nmatrix& init_val, ntype threshold, eqmode eq_mode) const{
  return evolve_until_equilibrium_general(init_val, threshold, eq_mode);
}
