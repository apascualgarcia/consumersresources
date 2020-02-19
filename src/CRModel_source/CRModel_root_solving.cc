#include "CRModel.h"
#include <iostream>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_poly.h>
#include <array>

statistics compute_critical_Delta(Metaparameters metaparams, stabilitymode stab_mode){
  delta_solver solv_params = {fitmode(sigmoidal), metaparams.equilibrium, stab_mode};
  return compute_critical_Delta(metaparams, solv_params);
}

nvector find_rough_interval_polynomial_fit(gsl_function* f, unsigned int Npoints, unsigned int verbose, interval bounds){
  nvector interval(Npoints, 0.);
  try{
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver* s;

    double x_lo = bounds.begin, x_hi =bounds.end, r=0.1;
    int status;
    int iter = 0;

    if(GSL_FN_EVAL(f, x_lo)*GSL_FN_EVAL(f,x_hi)>0){
      error err("It seems like the bounds of the initial interval in find_rough_interval_polynomial_fit do not contain zero. Return interval made of zeros.",1);
      throw err;
    }

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);

    gsl_root_fsolver_set(s, f, x_lo, x_hi);

    // the idea is first to find an interval where the solution roughly should be
    do{
      iter++;
      status = gsl_root_fsolver_iterate(s);
      r = gsl_root_fsolver_root(s);
      x_lo = gsl_root_fsolver_x_lower(s);
      x_hi = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(x_lo, x_hi, 0, 0.1);
      if(status==GSL_SUCCESS and verbose > 0){
        std::cout << "Found an interval for Delta critical : [" << x_lo << ";" << x_hi <<"]" << std::endl;
      }
    }while(status==GSL_CONTINUE);

    gsl_root_fsolver_free(s);

    size_t interval_length = Npoints;

    for(size_t i=0; i < interval_length; ++i){
      interval[i]=(x_lo+i*(x_hi-x_lo)/(interval_length-1));
    }
  }catch(error e){
    e.handle();
  }
  return interval;
}
nvector find_rough_interval_sigmoidal_fit(gsl_function* f, unsigned int Npoints, unsigned int verbose, interval bounds){
  nvector rough_interval;
  Solver_Parameters* s = (Solver_Parameters*) f->params;
  double initial_target = s->target;

  /* the idea is to find the first point x0 such that f(x0) > 0. It will be taken as the start of the interval */
  if(verbose > 0){
    std::cout << "We look for a point slightly above zero." << std::endl;
  }
  s->target = 0.01;
  double x_lo = find_zero(f, verbose, bounds);
  /* we do not forget to set the target back to its initial value */
  s->target = initial_target;


  /* we then find a second point such that f(x1) < 1. It will be taken as the end of the interval */
  if(verbose > 0){
    std::cout << "We then look for a high point slightly below one." << std::endl;
  }

  bool high_point_greater_than_one = false;
  /*  we first check if the point where f(x1) = 1-eps is outside of [0,1].
      This happens if f(1) < 1., we check that here */
  double x_hi=1.;
  s->target=0.;
  double proba_extinction_at_one = GSL_FN_EVAL(f,1.);
  s->target=initial_target;
  if(proba_extinction_at_one < 1.){
    std::cout << " Probability of extinction at one is smaller than one : " << proba_extinction_at_one<< std::endl;
    high_point_greater_than_one = true;
    if(verbose>0){
      std::cout << "It turns out that one does not saturate, so we take it as the end of the interval" << std::endl;
    }
  }


  if(not(high_point_greater_than_one)){
    s->target = 0.9;
    x_hi = find_zero(f, verbose, bounds);
    /* we do not forget to set the target back to its initial value */
    s->target = initial_target;
  }



  interval r_interval(x_lo, x_hi);
  if(verbose > 0){
    std::cout << "Found an interval for the critical value :" << r_interval<< std::endl;
  }

  for(size_t i = 0; i < Npoints; ++i){
    rough_interval.push_back(r_interval.begin + i*(r_interval.end-r_interval.begin)/(Npoints-1.));
  }
  return rough_interval;
}

double find_zero(gsl_function* f, unsigned int verbose, interval bounds){
  double estimate = 0.;
  try{
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver* s;

    double x_lo = bounds.begin, x_hi = bounds.end, r=0.05;
    double tolerance = 0.01;
    int status;
    int iter = 0;

    if(GSL_SIGN(GSL_FN_EVAL(f, x_lo)*GSL_FN_EVAL(f, x_hi))>0){
      estimate=NUMERICAL_ERROR;
      error err("It looks like the function does not vanish in the promised interval. Will return the error value",1);
      throw err;
    }

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);

    gsl_root_fsolver_set(s, f, x_lo, x_hi);

    do{
      iter++;
      status = gsl_root_fsolver_iterate(s);
      r = gsl_root_fsolver_root(s);
      x_lo = gsl_root_fsolver_x_lower(s);
      x_hi = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(x_lo, x_hi, 0, tolerance);
    }while(status==GSL_CONTINUE);
    estimate = r;
    if((estimate<x_lo) || (estimate>x_hi)){
      std::cout << "estimate outside of bound" << std::endl;
      estimate=NUMERICAL_ERROR;
      error err("It looks like the proposed zero is not in the desired interval. Will return the error value",1);
      throw err;
    }
    gsl_root_fsolver_free(s);
  }catch(error e){
    e.handle();
  }
  return estimate;
}
nvector find_rough_interval(gsl_function* f, unsigned int Npoints, unsigned int verbose, fitmode fit_mode, interval bounds){
  switch(fit_mode){
    case polynomial:
      return find_rough_interval_polynomial_fit(f, Npoints, verbose, bounds);
      break;

    case sigmoidal:
      return find_rough_interval_sigmoidal_fit(f, Npoints, verbose, bounds);
      break;

    case sigmoidal_erf:
      return find_rough_interval_sigmoidal_fit(f, Npoints, verbose, bounds);
      break;

    default:
      error err("Cannot find rought interval to compute delta critical with that fitting mode. It either does not exist or has not been implemented yet.");
      throw err;
      break;
  }
}
statistics compute_critical_Delta(Metaparameters metaparams, delta_solver delta_solv){
  eqmode equilibrium = delta_solv.eq_mode;
  unsigned int Nsimul_frun;
  unsigned int Nsimul_srun;
  size_t interval_length;
  double target;
  interval initial_guess(0., 1.);

  /* tells what is the target for the function we want to solve */
  switch(delta_solv.eq_mode){
    case oneextinct:
    /* in the case 'one extinction' we want the probability to be equal to 0.5 */
      target=0.5;
      break;
    case convergence:
    /* in the case 'convergence' we want the average number of extinctions to be 1 */
      target=1.;
      break;
    default:
      error err("This type of equilibrium has not been implemented yet or does not exist.");
      throw err;
      break;
  }

  switch(delta_solv.fit_mode){
    case polynomial:
      Nsimul_frun = COMPUTE_CRITICAL_DELTA_POLY_FIRST_RUN;
      Nsimul_srun = COMPUTE_CRITICAL_DELTA_POLY_SECOND_RUN;
      interval_length = COMPUTE_CRITICAL_DELTA_POLY_INTERVAL_LENGTH;
      break;
    case sigmoidal:
      Nsimul_frun = COMPUTE_CRITICAL_DELTA_SIG_FIRST_RUN;
      Nsimul_srun = COMPUTE_CRITICAL_DELTA_SIG_SECOND_RUN;
      interval_length= COMPUTE_CRITICAL_DELTA_SIG_INTERVAL_LENGTH;
      break;
    default:
      error err("This type of fitting mode has not been implemented yet or does not exist.");
      throw err;
      break;
  }


  Solver_Parameters params;
  params.metaparameters = &metaparams;
  params.Nsimul = Nsimul_frun;
  params.equilibrium = equilibrium;
  params.target = target;

  gsl_function F;
  F.function = &function_av_extinct_solver;
  F.params = &params;

  if(metaparams.verbose > 0){
    std::cout << "Now attempting to find the critical delta for the following set of parameters : " << metaparams << std::endl;
    std::cout << "We first find a rough interval where we know the critical delta will lie ("<<params.Nsimul << " runs per point). " << std::endl;
    std::cout << "We first test if the probability of extinction is strictly greater than zero for Delta=1" << std::endl;
  }

  double proba_extinction_at_one = GSL_FN_EVAL(&F, 1.);
  if(proba_extinction_at_one<=0){
    if(metaparams.verbose > 0){
      std::cout << "We will hence return the default value of Delta=1 which is synonym of global stability." << std::endl;
    }
    statistics globally_stable;
    globally_stable.mean_=1.;
    return globally_stable;
  }else if(proba_extinction_at_one<0.5){
    if(metaparams.verbose > 0){
      std::cout << "We expect the critical delta to be outside the range [0,1] " << std::endl;
    }
  }
  nvector interval = find_rough_interval(&F, interval_length, metaparams.verbose, delta_solv.fit_mode, initial_guess);

  /* when we have found an interval we highly suspect of containing the root
    we compute the average number of extinctions at a better accuracy for the points
    in the interval */
  params.Nsimul=Nsimul_srun;

  nvector function_y_values;
  if(metaparams.verbose > 0){
    std::cout << "Now computing the ";
    if(equilibrium == convergence){
      std::cout << "average number of extinctions";
    }else if(equilibrium==oneextinct){
      std::cout << "probability of getting more than one extinction";
    }
    std::cout << " for "<< interval_length << " points inside this interval (" << params.Nsimul <<" runs per point)" << std::endl;
  }

  for(size_t i=0; i < interval_length; ++i){
    double result = function_av_extinct_solver(interval[i], &params);
    function_y_values.push_back(result);
  }

  /* after getting these ten points, we fit them with a curve of a given shape,
   that allows us to estimate delta critical */
  statistics delta_crit = estimate_delta_crit_from_interval(interval, function_y_values, metaparams, delta_solv);
  return delta_crit;
}
statistics compute_critical_alpha(Metaparameters& metaparams, ntype accuracy, fitmode fit_mode){
  statistics critical_alpha;
  double target=0.;
  unsigned int Nsimul_frun = 100, Nsimul_srun = 10000;
  size_t interval_length = 50;
  interval initial_guess(0., metaparams.physical_maximum_alpha0());

  Solver_Parameters params;
  params.metaparameters = &metaparams;
  params.Nsimul = Nsimul_frun;
  params.target = target;

  gsl_function F;
  F.function = &function_proba_feasability_solver;
  F.params = &params;

  if(metaparams.verbose > 0){
    std::cout << "Now attempting to find the critical alpha for the following set of parameters : " << metaparams << std::endl;
    std::cout << "We first find a rough interval where we know the critical alpha will lie ("<<params.Nsimul << " runs per point). " << std::endl;
  }
  nvector interval = find_rough_interval(&F, interval_length, metaparams.verbose, fit_mode, initial_guess);
  /* when we have found an interval we highly suspect of containing the root
    we compute the feasability probability at a better accuracy for the points
    in the interval */
  params.Nsimul=Nsimul_srun;

  nvector function_y_values;
  if(metaparams.verbose > 0){
    std::cout << "Now computing the feasability probability for different alpha, ";
    std::cout << interval_length << " points inside this interval (" << params.Nsimul <<" runs per point)" << std::endl;
  }

  for(size_t i=0; i < interval_length; ++i){
    double result = function_proba_feasability_solver(interval[i], &params);
    function_y_values.push_back(result);
  }

  /* after getting these ten points, we fit them with a curve of a given shape,
   that allows us to estimate the critical alpha */
  critical_alpha = estimate_alpha_crit_from_interval(interval, function_y_values, metaparams, fit_mode);
  return critical_alpha;
}


double function_av_extinct_solver(double delta, void*params){
  Solver_Parameters* s = (Solver_Parameters*) params;
  Metaparameters* m = s->metaparameters;
  unsigned int Nsimul = s->Nsimul;
  double target = s->target;

  switch(s->equilibrium){
    case convergence:
      return average_number_of_extinctions(delta, m, Nsimul)-target;
      break;
    case oneextinct:
      return probability_of_extinction_greather_than_one(m, delta, Nsimul)-target;
      break;
    default:
      error err("Equilibrium type not implemented yet.");
      throw err;
      break;
  }
}
double function_proba_feasability_solver(double alpha, void* params){
  Solver_Parameters* s = (Solver_Parameters*) params;
  Metaparameters* m = s->metaparameters;
  m->alpha0 = alpha;
  unsigned int Nsimul = s->Nsimul;
  double target = s->target;

  return find_feasability_probability(*m, Nsimul)-target;
}

double function_proba_feasability_solver_S0(double S, void* params){
  Solver_Parameters* s = (Solver_Parameters*) params;
  Metaparameters* m = s->metaparameters;
  m->S0 = S;
  unsigned int Nsimul = s->Nsimul;
  double target = s->target;

  return find_feasability_probability(*m, Nsimul)-target;

}

statistics solve_for_delta_with_fit(fitting_parameters& fit_parameters, double & x_lo, double & x_hi, const Metaparameters& m, delta_solver delta_solv){
  double estimate = 0., error = 0.;
  if(m.verbose > 0){
    std::cout << "Now we find the zero of the fit to determine delta critical (parameters =";
    for(size_t i = 0; i < fit_parameters.fit_parameters->size; ++i){
      std::cout << " "<<gsl_vector_get(fit_parameters.fit_parameters, i) << "+/-";
      std::cout << gsl_vector_get(fit_parameters.error, i);
    }
    std::cout << ")" << std::endl;
  }

  switch(delta_solv.fit_mode){
    case sigmoidal:{
      estimate = gsl_vector_get(fit_parameters.fit_parameters,0);
      error = gsl_vector_get(fit_parameters.error, 0);
      break;
    }
    case sigmoidal_erf:{
      estimate = gsl_vector_get(fit_parameters.fit_parameters,0);
      error = gsl_vector_get(fit_parameters.error, 0);
    }
    case polynomial:{
      struct error err("PLEASE IMPLEMENT THE FUNCTION TO SOLVE FOR DELTA WITH POLYNOMIAL FIT.");
      throw err;
      break;
    }

    default:{
      struct error err("This type of fitting mode has not been implemented yet or does not exist.");
      throw err;
      break;
    }
  }
  statistics delta;
  delta.mean_ = estimate;
  delta.std_deviation_ = error;
  delta.median_ = NULL;
  return delta;
}
