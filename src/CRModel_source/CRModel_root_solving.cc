#include "CRModel.h"
#include <iostream>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_poly.h>
#include <array>

statistics compute_critical_Delta(Metaparameters metaparams, ntype accuracy){
  delta_solver solv_params = {fitmode(sigmoidal), metaparams.equilibrium};
  return compute_critical_Delta(metaparams, accuracy, solv_params);
}
nvector find_rough_interval_polynomial_fit(gsl_function* f, unsigned int Npoints, unsigned int verbose){
  nvector interval;

  const gsl_root_fsolver_type *T;
  gsl_root_fsolver* s;

  double x_lo = 0., x_hi =1., r=0.1;
  int status;
  int iter = 0;

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
    interval.push_back(x_lo+i*(x_hi-x_lo)/(interval_length-1));
  }
  return interval;
}
nvector find_rough_interval_sigmoidal_fit(gsl_function* f, unsigned int Npoints, unsigned int verbose){
  nvector interval;
  Solver_Parameters* s = (Solver_Parameters*) f->params;
  double initial_target = s->target;

  /* the idea is to find the first point x0 such that f(x0) > 0. It will be taken as the start of the interval */
  if(verbose > 0){
    std::cout << "We look for a low point such that we have a non-zero chance of observing an extinction." << std::endl;
  }
  s->target = 0.01;
  double x_lo = find_zero(f, Npoints, verbose);

  /* we then find a second point such that f(x1) < 1. It will be taken as the end of the interval */
  if(verbose > 0){
    std::cout << "We then look for a high point such that we will not always observe at least an extinction" << std::endl;
  }
  s->target = 0.8;
  double x_hi = find_zero(f, Npoints, verbose);

  /* we do not forget to set the target back to its initial value */
  s->target = initial_target;

  if(verbose > 0){
    std::cout << "Found an interval for Delta critical : [" << x_lo << ";" << x_hi <<"]" << std::endl;
  }

  for(size_t i = 0; i < Npoints; ++i){
    interval.push_back(x_lo + i*(x_hi-x_lo)/(Npoints-1.));
  }
  return interval;
}

double find_zero(gsl_function* f, unsigned int Npoints, unsigned int verbose){
  double estimate = 0.;

  const gsl_root_fsolver_type *T;
  gsl_root_fsolver* s;

  double x_lo = 0., x_hi =1., r=0.05;
  double tolerance = 0.05;
  int status;
  int iter = 0;

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
    status = gsl_root_test_interval(x_lo, x_hi, 0, tolerance);
  }while(status==GSL_CONTINUE);
  estimate = r;

  gsl_root_fsolver_free(s);
  return estimate;
}

nvector find_rough_interval(gsl_function* f, unsigned int Npoints, unsigned int verbose, fitmode fit_mode){
  switch(fit_mode){
    case polynomial:
      return find_rough_interval_polynomial_fit(f, Npoints, verbose);
      break;

    case sigmoidal:
      return find_rough_interval_sigmoidal_fit(f, Npoints, verbose);
      break;

    default:
      std::cerr << "Cannot find rought interval to compute delta critical with that fitting mode" << std::endl;
      std::cerr << "It either does not exist or has not been implemented yet. Aborting the simulation now" << std::endl;
      abort();
  }
}

statistics compute_critical_Delta(Metaparameters metaparams, ntype accuracy, delta_solver delta_solv){
  eqmode equilibrium = delta_solv.eq_mode;
  unsigned int Nsimul_frun;
  unsigned int Nsimul_srun;
  size_t interval_length;
  double target;

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
      std::cerr << "This type of equilibrium has not been implemented yet or does not exist " << std::endl;
      std::cerr << "Abort the program now " << std::endl;
      abort();
      break;
  }

  switch(delta_solv.fit_mode){
    case polynomial:
      Nsimul_frun = 50;
      Nsimul_srun = 100;
      interval_length = 10;
      break;
    case sigmoidal:
      Nsimul_frun = 50;
      Nsimul_srun = 50;
      interval_length= 50;
      break;
    default:
      std::cerr << "This type of fitting mode has not been implemented yet or does not exist" << std::endl;
      std::cerr << "Aborting the program now" << std::endl;
      abort();
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
  }

  nvector interval = find_rough_interval(&F, interval_length, metaparams.verbose, delta_solv.fit_mode);

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
      std::cerr << "equilibrium type not implemented yet"<< std::endl;
      std::cerr << "Aborting simulation now" << std::endl;
      abort();
      break;
  }

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
    case sigmoidal:
      estimate = gsl_vector_get(fit_parameters.fit_parameters,0);
      error = gsl_vector_get(fit_parameters.error, 0);
      break;
    case polynomial:
      std::cerr << "PLEASE IMPLEMENT THE FUNCTION TO SOLVE FOR DELTA WITH POLYNOMIAL FIT" << std::endl;
      std::cerr << "Aborting the simulation now";
      abort();
      break;

    default:
      std::cerr << "This type of fitting mode has not been implemented yet or does not exist" << std::endl;
      std::cerr << "Abort simulation now " << std::endl;
      abort();
      break;
  }
  statistics delta = {estimate, error, NULL};
  return delta;
}
