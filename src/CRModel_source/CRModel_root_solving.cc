#include "CRModel.h"
#include <iostream>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <array>

Delta_critical compute_critical_Delta(Metaparameters metaparams, ntype accuracy){
  Delta_critical delta_crit;

  // set up the solver
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver* s;

  double x_lo = 0., x_hi =1., r=0.1;
  int status;
  int iter = 0;

  Solver_Parameters params;
  params.metaparameters = &metaparams;
  params.Nsimul = 10;

  gsl_function F;
  F.function = &function_av_extinct_solver;
  F.params = &params;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set(s, &F, x_lo, x_hi);

  // The idea is to find a smaller and smaller interval in which we have a higher and higher accuracy
  // i.e. you narrow down the interval
  unsigned int Nsimul_max = 1./accuracy;
  std::array<unsigned int,2> Nsimul_for_intervals={100, Nsimul_max};
  std::array<ntype,2> accuracy_for_intervals={0.1, accuracy};

  for(size_t i = 0; i < Nsimul_for_intervals.size(); ++i){
    params.Nsimul = Nsimul_for_intervals[i];
    F.params = &params;
    gsl_root_fsolver_set(s, &F, x_lo, x_hi);
    do{
      iter++;
      status = gsl_root_fsolver_iterate(s);
      r = gsl_root_fsolver_root(s);
      x_lo = gsl_root_fsolver_x_lower(s);
      x_hi = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(x_lo, x_hi, 0, accuracy_for_intervals[i]);
      if(status==GSL_SUCCESS){
        std::cout << "Found an interval for Delta critical : [" << x_lo << ";" << x_hi <<"]" << std::endl;
        std::cout << "In this interval we have Dcrit=" << r << " with number of extinctions " << average_number_of_extinctions(r, &params) << std::endl;
      }
    }while(status==GSL_CONTINUE);
    // when we have found an interval good enough for the desired accuracy, we refine it and increase the accuracy
  }

  gsl_root_fsolver_free(s);

  delta_crit.delta_crit = r;
  delta_crit.delta_low = x_lo;
  delta_crit.delta_high = x_hi;
  delta_crit.accuracy = accuracy;

  return delta_crit;
}


double average_number_of_extinctions(double delta, void* params){
  Solver_Parameters* s = (Solver_Parameters*) params;
  Metaparameters* m = s->metaparameters;
  unsigned int Nsimul = s->Nsimul;

  Extinction ext = compute_average_extinction(m, ntype(delta), Nsimul);
  double av_number_extinct = double(ext.extinct);

  return av_number_extinct;
}

double function_av_extinct_solver(double delta, void*params){
  return average_number_of_extinctions(delta, params)-1.;
}
