#include "CRModel.h"
#include <iostream>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_poly.h>
#include <array>

double compute_critical_Delta(Metaparameters metaparams, ntype accuracy){
  return compute_critical_Delta(metaparams, accuracy, metaparams.equilibrium);
}

nvector find_rough_interval(gsl_function* f, unsigned int Npoints, unsigned int verbose){
  nvector interval;

  const gsl_root_fsolver_type *T;
  gsl_root_fsolver* s;

  double x_lo = 0., x_hi =1., r=0.1;
  int status;
  int iter = 0;

  gsl_function F = *f;

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

double compute_critical_Delta(Metaparameters metaparams, ntype accuracy, eqmode equilibrium){
  double delta_crit=0.;
  unsigned int Nsimul_frun = 100;
  unsigned int Nsimul_srun = 1000;
  size_t interval_length = 10;

  Solver_Parameters params;
  params.metaparameters = &metaparams;
  params.Nsimul = Nsimul_frun;
  params.equilibrium = equilibrium;

  gsl_function F;
  F.function = &function_av_extinct_solver;
  F.params = &params;

  if(metaparams.verbose > 0){
    std::cout << "Now attempting to find the critical delta for the following set of parameters : " << metaparams << std::endl;
    std::cout << "We first find a rough interval where we know the critical delta will lie. " << std::endl;
  }

  nvector interval = find_rough_interval(&F, interval_length, metaparams.verbose);

  if(metaparams.verbose > 0){
    std::cout << "Now computing the ";
    if(equilibrium == convergence){
      std::cout << "average number of extinctions";
    }else if(equilibrium==oneextinct){
      std::cout << "probability of getting more than one extinction";
    }
    std::cout << " for ten points inside this interval (" << params.Nsimul <<" runs per point)" << std::endl;
  }

  /* when we have found an interval we highly suspect of containing the root
    we compute the average number of extinctions at a better accuracy for the points
    in the interval */
  params.Nsimul=Nsimul_srun;
  F.params = &params;
  nvector function_y_values;

  for(size_t i=0; i < interval_length; ++i){
    double result = function_av_extinct_solver(interval[i], &params);
    function_y_values.push_back(result);
  }

  /* after getting these ten points, we fit them with a curve of a given shape,
   that allows us to estimate delta critical */
  delta_crit = estimate_delta_crit_from_interval(interval, function_y_values, metaparams, equilibrium);

  return delta_crit;
}




double function_av_extinct_solver(double delta, void*params){
  Solver_Parameters* s = (Solver_Parameters*) params;
  Metaparameters* m = s->metaparameters;
  unsigned int Nsimul = s->Nsimul;

  if(s->equilibrium==convergence){
    return average_number_of_extinctions(delta, m, Nsimul)-1.;
  }else if(s->equilibrium==oneextinct){
    return probability_of_extinction_greather_than_one(m, delta, Nsimul)-0.5;
  }else{
    std::cerr << "equilibrium type not implemented yet, return 0"<< std::endl;
    return 0.;
  }
}

double solve_for_delta_with_fit(const gsl_vector* fit_parameters, double & x_lo, double & x_hi, const Metaparameters& m, eqmode equilibrium){
  double estimate = 0.;

  if(m.verbose > 0){
    std::cout << "Now we find the zero of the fit to determine delta critical (parameters =";
    for(size_t i = 0; i < fit_parameters->size; ++i){
      std::cout << " "<<gsl_vector_get(fit_parameters, i);
    }
    std::cout << ")" << std::endl;
  }

  /*
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
  std::cout << "All good until here" << std::endl;
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
  */

  if(m.verbose > 0){
    std::cout << "We interpret the fitting parameters as coefficients of a degree 3 polynomial (which should be true, please check that)" << std::endl;
  }

  /* initiated these as -10 because we know it can't be a root */
  double x0 = -10;
  double x1 = -10;
  double x2 = -10;

  double a0 = gsl_vector_get(fit_parameters,0), a1= gsl_vector_get(fit_parameters,1), a2 = gsl_vector_get(fit_parameters,2),a3 = gsl_vector_get(fit_parameters,3);
  double c = a0/a3, b = a1/a3, a = a2/a3;

  int success = gsl_poly_solve_cubic(a,b,c, &x0, &x1, &x2);
  if(x0 == -10){
    std::cerr << "Error in the root finding, could not find a single root" << std::endl;
  }
  if(x1==-10){
    estimate = x0;
  }else{
    unsigned int Nsimul =100;
    if(m.verbose > 0){
      std::cout << "Found three potential roots, compute which one is the best (" << Nsimul <<" runs per point)" << std::endl;
      std::cout << "Please implement the different equilibrium modes in case of three roots" << std::endl;
    }
    estimate = -1.;
    /*
    Metaparameters m_copy = m;
    Extinction_statistics av0 = compute_average_extinction(&m_copy, ntype(x0), Nsimul);
    Extinction_statistics av1 = compute_average_extinction(&m_copy, ntype(x1), Nsimul);
    Extinction_statistics av2 = compute_average_extinction(&m_copy, ntype(x2), Nsimul);

    double dist0 = (av0.extinct.mean-1.)*(av0.extinct.mean-1.);
    double dist1 = (av1.extinct.mean-1.)*(av1.extinct.mean-1.);
    double dist2 = (av2.extinct.mean-1.)*(av2.extinct.mean-1.);

    if(dist0 < dist1){
      if(dist0 < dist2){
        estimate = x0;
      }else{
        estimate = x2;
      }
    }else{
      if(dist1 < dist2){
        estimate = x1;
      }else{
        estimate = x2;
      }
    }
    */
  }
  return estimate;
}
