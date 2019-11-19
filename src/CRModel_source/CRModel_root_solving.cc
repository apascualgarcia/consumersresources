#include "CRModel.h"
#include <iostream>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <array>

double compute_critical_Delta(Metaparameters metaparams, ntype accuracy){
  double delta_crit=0.;

  // set up the solver
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver* s;

  double x_lo = 0., x_hi =1., r=0.1;
  int status;
  int iter = 0;

  Solver_Parameters params;
  params.metaparameters = &metaparams;
  params.Nsimul = 50;

  gsl_function F;
  F.function = &function_av_extinct_solver;
  F.params = &params;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);

  gsl_root_fsolver_set(s, &F, x_lo, x_hi);

  // the idea is first to find an interval where the solution roughly should be
  do{
    iter++;
    status = gsl_root_fsolver_iterate(s);
    r = gsl_root_fsolver_root(s);
    x_lo = gsl_root_fsolver_x_lower(s);
    x_hi = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(x_lo, x_hi, 0, 0.1);
    if(status==GSL_SUCCESS and metaparams.verbose > 0){
      std::cout << "Found an interval for Delta critical : [" << x_lo << ";" << x_hi <<"]" << std::endl;
    }
  }while(status==GSL_CONTINUE);

  /* when we have found an interval we highly suspect of containing the root
    we compute the average number of extinctions at a better accuracy for ten points
   in the interval */
  nvector interval;
  size_t interval_length = 10;
  for(size_t i=0; i < interval_length; ++i){
    interval.push_back(x_lo+i*(x_hi-x_lo)/(interval_length-1));
  }

  params.Nsimul=500;
  F.params = &params;
  nvector extinctions;
  for(size_t i=0; i < interval_length; ++i){
    double result = function_av_extinct_solver(interval[i], &params);
    extinctions.push_back(result);
  }

  /* after getting these ten points, we fit them with a curve of a given shape,
   that allows us to estimate delta critical */
  delta_crit = estimate_delta_crit_from_interval(interval, extinctions);
  gsl_root_fsolver_free(s);

  return delta_crit;
}


double average_number_of_extinctions(double delta, void* params){
  Solver_Parameters* s = (Solver_Parameters*) params;
  Metaparameters* m = s->metaparameters;
  unsigned int Nsimul = s->Nsimul;

  Extinction_statistics ext = compute_average_extinction(m, ntype(delta), Nsimul);
  double av_number_extinct = double(ext.extinct.mean);

  return av_number_extinct;
}

double function_av_extinct_solver(double delta, void*params){
  return average_number_of_extinctions(delta, params)-1.;
}

double solve_for_delta_with_fit(const gsl_vector* fit_parameters, double & x_lo, double & x_hi){
  double estimate = 0.;

  unsigned int max_iter = 100;
  unsigned int N = NUMBER_OF_FITTING_PARAMETERS;
  double params[N];
  for(size_t i = 0; i < N; ++i){
    params[i] = gsl_vector_get(fit_parameters,i);
  }

  int status;
  int iter = 0;

  const gsl_root_fdfsolver_type* T;
  gsl_root_fdfsolver* s;
  double r = 0.5*(x_lo+x_hi);
  double x0;
  gsl_function_fdf F;

  F.f = &choice_of_fitting_function;
  F.df = NULL;
  F.fdf = NULL;
  F.params = &params;

  T = gsl_root_fdfsolver_newton;
  s = gsl_root_fdfsolver_alloc(T);
  gsl_root_fdfsolver_set(s, &F, r);
  do{
    iter++;
    status = gsl_root_fdfsolver_iterate(s);
    x0 = r;
    status = gsl_root_test_delta(r, x0, 0, 1e-3);

  }while(status==GSL_CONTINUE && iter < max_iter);
  if(iter > max_iter){
    std::cerr << "The solver wasn't able to estimate delta critical, the best estimate will be returned" << std::endl;
  }
  estimate = r;
  gsl_root_fdfsolver_free(s);
  return estimate;
}
