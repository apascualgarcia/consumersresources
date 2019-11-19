#include "CRModel.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <iostream>
#include <array>
/*
    code from this file is heavily inspired from the GNU docs example of function fitting
    https://www.gnu.org/software/gsl/doc/html/nls.html#examples
*/

struct data {
  size_t p;
  size_t n;
  double* x;
  double* y;
};

/*  That's the choice of the function we use for fitting, i.e. assume Yi = f(Xi)
    We choose Yi = a0+a1*Xi+a2*Xi^2+a3*Xi^3 as the function to fit  */
double choice_of_fitting_function(double x, void* params){
  double* a = (double*)params;
  double result = 0.;
  unsigned int N = NUMBER_OF_FITTING_PARAMETERS;
  for(size_t i = 0; i < N; ++i){
    result+= a[i]*gsl_pow_uint(x,i);
  }
  return result;
}

double fitting_function(double x, const gsl_vector* a){
  unsigned int N = NUMBER_OF_FITTING_PARAMETERS;
  double params[N];
  for(size_t i=0; i < N; ++i){
    params[i] = gsl_vector_get(a, i);
  }
  return choice_of_fitting_function(x, &params);
}


int function_to_fit(const gsl_vector*  params, void* data, gsl_vector* f){
  size_t n = ((struct data*)data)->n;
  double* x = ((struct data*)data)->x;
  double* y = ((struct data*)data)->y;

  for(size_t i=0; i < n; ++i){
    double Yi = fitting_function(x[i], params);
    gsl_vector_set(f,i,Yi-y[i]);
  }

  return GSL_SUCCESS;
}
void callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w){
  gsl_vector *f = gsl_multifit_nlinear_residual(w);
  gsl_vector *x = gsl_multifit_nlinear_position(w);
  double rcond;

  /* compute reciprocal condition number of J(x) */
  gsl_multifit_nlinear_rcond(&rcond, w);

  fprintf(stderr, "iter %2zu: a3 = %.4f, a2 = %.4f, a1 = %.4f, cond(J) = %8.4f, |f(x)| = %.4f\n",
          iter,
          gsl_vector_get(x, 0),
          gsl_vector_get(x, 1),
          gsl_vector_get(x, 2),
          1.0 / rcond,
          gsl_blas_dnrm2(f));
  return;
}
void fit_points_with_function(const nvector& interval, const nvector& points, gsl_vector* fit_parameters){

  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_workspace *w;
  gsl_multifit_nlinear_fdf fdf;
  gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();

  const unsigned int N = interval.size();

  const size_t n = N;
  const size_t p = NUMBER_OF_FITTING_PARAMETERS;

  gsl_vector* f;
  double x[N], y[N], weights[N];
  struct data d = {p, n, x, y};
  double a_init[p];

  for(size_t i=0; i < p; ++i){
    a_init[i]=1.0;
  }

  gsl_vector_view a = gsl_vector_view_array(a_init, p);
  gsl_vector_view wts = gsl_vector_view_array(weights, n);
  gsl_rng* r;
  int status, info;
  size_t i;

  const double atol = 1e-8;
  const double gtol = 1e-8;
  const double ftol = 0.0;

  gsl_rng_env_setup();
  r = gsl_rng_alloc(gsl_rng_default);


  /* define the function to be minimized*/
  fdf.f = function_to_fit;
  fdf.df = NULL;
  fdf.fvv = NULL;
  fdf.n = n;
  fdf.p = p;
  fdf.params = &d;

  /* this is the data to be fitted */
  for (i=0; i < n; ++i){
    x[i] = double(interval[i]);
    y[i] = double(points[i]);
    weights[i] = 1.;
  }

  /* allocate workspace with default parameters */
  w = gsl_multifit_nlinear_alloc(T, &fdf_params, n, p);

  /* initialize solver with starting points and weights */
  gsl_multifit_nlinear_winit(&a.vector, &wts.vector, &fdf, w);

  /* compute initial cost function */
  double chisq0;
  f = gsl_multifit_nlinear_residual(w);
  gsl_blas_ddot(f,f, &chisq0);

  /* solve the system with a maximum of 100 iterations */
  status = gsl_multifit_nlinear_driver(100, atol, gtol, ftol, NULL, NULL, &info, w);

  for(size_t i=0; i < w->x->size; ++i){
    gsl_vector_set(fit_parameters, i, gsl_vector_get(w->x,i));
  }

  gsl_multifit_nlinear_free(w);
  gsl_rng_free(r);

  return;
}

double estimate_delta_crit_from_interval(const nvector& interval, const nvector& extinctions, const Metaparameters& m){
  double delta_crit=0.;
  double x_lo = interval[0];
  double x_hi = interval[interval.size()-1];

  /*so we have this interval of x-points with their y-values, we choose to make a
   polynomial fit around it and take delta critical as the root of that polynomial */

  /* first select the points to fit */
  nvector x_points_to_fit = interval;
  nvector y_points_to_fit = extinctions;

  /* Then we actually find the parameters that fit our choice of function best */
  unsigned int number_of_fitting_parameters = NUMBER_OF_FITTING_PARAMETERS;
  gsl_vector* fit_parameters = gsl_vector_alloc(number_of_fitting_parameters);
  if(m.verbose > 0){
    std::cout << "Now fitting the " << x_points_to_fit.size() << " chosen into the specific function" << std::endl;
  }
  fit_points_with_function(x_points_to_fit, y_points_to_fit, fit_parameters);


  /* with the fitting parameters estimated, we can actually solve for Delta numerically */
  delta_crit = solve_for_delta_with_fit(fit_parameters, x_lo, x_hi, m);
  gsl_vector_free(fit_parameters);
  std::cout << "zero estimated at " << delta_crit << std::endl;

  /* finally, we return the estimated value */
  return delta_crit;
}
