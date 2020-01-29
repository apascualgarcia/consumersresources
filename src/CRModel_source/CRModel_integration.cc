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


Extinction CRModel::evolve_until_equilibrium_general(const nmatrix& init_val, ntype threshold, eqmode eq_mode, writemode write_mode) const{
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
  gsl_odeiv2_system sys = {this->equations_of_evolution, NULL, dim, this->model_param->get_parameters()};

  unsigned int counts = 0;
  double previous_y[p->NR+p->NS][INDICES_FOR_AVERAGE];
  double diff[p->NR+p->NS][INDICES_FOR_AVERAGE];
  bool numerical_convergence = false;
  bool exit_condition = false;
  bool one_extinct = false;
  /* if it turns out that the system was ALREADY at equilibrium, we will return 0 as the final time */
  bool started_at_equilibrium = false;
  std::vector<unsigned int> extinct_variables;
  double t_eq = 0.;

  exit_condition = (t>=tmax) or numerical_convergence or started_at_equilibrium;

  /* a rough estimate (higher than the actual value) on the error of the solution*/


  if(write_mode.write){
    write_mode.write_path << std::setprecision(print_precision);
  }

  /* system temporal evolution */
  while(not(exit_condition)){

    /* store the previous ten values to estimate convergence */
    for(size_t i = 0; i < p->NR+p->NS; ++i){
      previous_y[i][counts%INDICES_FOR_AVERAGE] = y[i];
    }

    /* evolve the system to the next state */
    int status = gsl_odeiv2_evolve_apply(e,c,s, &sys, &t, tmax, &step_size, y);
    if (status != GSL_SUCCESS){
      std::cerr << "Error in the integration of the ODE!" << std::endl;
      abort();
      break;
    }
    ntype ERROR_ON_THE_SOLUTION =0.;
    for(size_t i=0; i < p->NR+p->NS;++i){
      ntype local_error=e->yerr[i];
      if(local_error > ERROR_ON_THE_SOLUTION){
        ERROR_ON_THE_SOLUTION = local_error;
      }
    }

    /* we write the system to the desired output if wanted */
    if(write_mode.write){
      write_mode.write_path << t << " ";
      for(size_t i =0; i<p->NR+p->NS; ++i){
        write_mode.write_path << y[i] << " ";
      }
      write_mode.write_path << ERROR_ON_THE_SOLUTION << " ";
    }

    /* we then compute the differences between this step and the previous one to check if we started already at equilibrium */
    if(counts <= INDICES_FOR_AVERAGE){
      for(size_t i=0; i < p->NR+p->NS; ++i){
        diff[i][counts%INDICES_FOR_AVERAGE] = y[i]-previous_y[i][counts%INDICES_FOR_AVERAGE];
      }
    }

    /* after ten counts we check whether the system was de facto at equilibrium */
    if(counts==INDICES_FOR_AVERAGE){
      started_at_equilibrium = true;
      for(size_t j=0; j < INDICES_FOR_AVERAGE and started_at_equilibrium; ++j){
        for(size_t k=0; (k < p->NR+p->NS) and started_at_equilibrium; ++k){
          double difference = abs(diff[k][j]);
          if(difference> INTEGRATOR_ZERO){
            started_at_equilibrium = false;
          }
        }
      }

      if(started_at_equilibrium and p->verbose > 3){
        std::cout << "\t \t \t It turns out the system was already at equilibrium!" << std::endl;
      }
    }

    /*  displays the values on cout if asked (this can also be done by writing the output to cout with
        writemode) */
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

    /* computes the "convergence coefficient" to estimate the convergence (one criterion to stop)*/
    numerical_convergence = convergence_criterion(threshold, t, previous_y,y,counts, this->model_param->get_parameters(), this->equations_of_evolution, write_mode);

    /* if a resource/consumer is too small, we effectively set it to zero */
    bool local_exit = false;
    for(size_t i=0; i < p->NR+p->NS and not(local_exit); ++i){
      if(y[i] < SPECIES_EXTINCT){
        y[i]=INTEGRATOR_ZERO;
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

    if(write_mode.write){
      write_mode.write_path << std::endl;
    }

    /* we check whether we should continue integrating or not */
    exit_condition = (t>=tmax) or numerical_convergence or started_at_equilibrium;
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
      }else if(numerical_convergence){
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

bool convergence_criterion(const double threshold, const double t, const double previous_y[][INDICES_FOR_AVERAGE], const double y[], unsigned int counts, void* params, func_equ_evol equ_evol, writemode wm){
  bool converged = false;
  Parameter_set* p = &(*(Parameter_set*) params);
  const unsigned int sys_size=p->NR+p->NS;
  if(counts>=INDICES_FOR_AVERAGE){
    double dydt[sys_size];
    /* computes dNi/dt in dydt */
    /* CHECK WHY THIS DOES NOT GIVE APPROPRIATE RESULTS */
    equ_evol(t, y, dydt, params);


    bool all_smaller_than_threshold=true;
    bool stay_condition = true;

    for(size_t i=0; stay_condition ;++i){
      double local_error=0.;
      /* we check the convergence only of non extinct species */
      if(y[i] > SPECIES_EXTINCT){
        local_error = abs(dydt[i]/y[i]);
      }
      if(local_error > threshold){
        all_smaller_than_threshold = false;
      }
      if(wm.write){
        wm.write_path << local_error << " ";
      }
      stay_condition=wm.write or all_smaller_than_threshold;
      /*  Add -1 because this is evaluated at the end of the loop so
          e.g. when sys_size=6, stay_condition should be wrong when i = 5*/
      stay_condition=stay_condition and i < sys_size-1;
    }
    converged = all_smaller_than_threshold;
  }else{
    if(wm.write){
      for(size_t i=0; i < sys_size;++i){
        wm.write_path << "NaN ";
      }
    }
  }

  return converged;
}

Extinction CRModel::evolve_until_equilibrium(ntype threshold, eqmode eq_mode, writemode write_mode) const{
  return evolve_until_equilibrium_general((*eq_vals)[0],threshold, eq_mode, write_mode);
}
Extinction CRModel::evolve_until_equilibrium_from_abundances(const nmatrix& init_val, ntype threshold, eqmode eq_mode, writemode write_mode) const{
  return evolve_until_equilibrium_general(init_val, threshold, eq_mode, write_mode);
}
