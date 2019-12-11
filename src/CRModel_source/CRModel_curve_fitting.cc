#include "CRModel.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_exp.h>
#include <iostream>
#include <array>
#include <string>
/*
    code from this file is heavily inspired from the GNU docs example of function fitting
    https://www.gnu.org/software/gsl/doc/html/nls.html#examples
*/

unsigned int NUMBER_OF_FITTING_PARAMETERS;

/*  That's the choice of the function we use for fitting, i.e. assume Yi = f(Xi)
    We choose Yi = a0+a1*Xi+a2*Xi^2+a3*Xi^3 as the function to fit  */

double polynomial_fit(double x, const gsl_vector* a){
  double result = 0.;
  for(size_t i = 0; i < NUMBER_OF_FITTING_PARAMETERS; ++i){
    result+= gsl_vector_get(a,i)*gsl_pow_uint(x,i);
  }
  return result;
}

double sigmoidal_fit(double x, const gsl_vector* a){
  double result = 0.;
  double x0 = gsl_vector_get(a,0);
  double scale = gsl_vector_get(a,1);

  /* CAREFUL : add -0.5 because that's the target */
  result = 1./(1.+gsl_sf_exp(-scale*(x-x0)))-0.5;
  return result;
}

int function_to_fit_df(const gsl_vector* params, void* data, gsl_matrix* J){
  set_jacobian_matrix_fit(params, data, J);
  return GSL_SUCCESS;
}

void set_jacobian_matrix_fit(const gsl_vector* params, void* data, gsl_matrix* J){
  fitmode fit_mode =((struct data*)data)->fit_mode;

  switch(fit_mode){
    case sigmoidal:
      set_jacobian_matrix_sigmoidal_fit(params, data, J);
      break;
    case polynomial:
      std::cerr << "Please implement the jacobian for the polynomial fit " << std::endl;
      std::cerr << "Aborting the simulation." << std::endl;
      abort();
      break;
    default:
      std::cerr << "The jacobian for this type of fitting has not been implemented or does not exist yet" << std::endl;
      std::cerr << "Aborting the simulation." << std::endl;
      abort();
      break;
  }
  return;
}

void set_jacobian_matrix_sigmoidal_fit(const gsl_vector* params, void* data, gsl_matrix* J){
  double* x = ((struct data*)data)->x;
  size_t n = ((struct data*)data)->n;

  double x0 = gsl_vector_get(params, 0);
  double scale = gsl_vector_get(params, 1);

  /* J has shape n x p i.e. here it has two columns for n rows */
  /* it is filled by df_i/dparams_j */
  for(size_t i = 0; i < n;++i){
    double factor = gsl_sf_exp(-scale*(x[i]-x0));
    factor /= gsl_pow_2(1.+gsl_sf_exp(-scale*(x[i]-x0)));
    gsl_matrix_set(J,i,0,-scale*factor);
    gsl_matrix_set(J,i,1,(x[i]-x0)*factor);
  }
  return;
}



double fitting_function(double x, const gsl_vector* a, fitmode fit_mode){
  switch(fit_mode){
    case sigmoidal:
      return sigmoidal_fit(x,a);
    case polynomial:
      return polynomial_fit(x,a);
    default:
      std::cerr << "This type of fitting either does not exist or has not been implemented yet " << std::endl;
      std::cerr << "Aborting the program now" << std::endl;
      abort();
      break;
  }
}

int function_to_fit(const gsl_vector*  params, void* data, gsl_vector* f){
  size_t n = ((struct data*)data)->n;
  double* x = ((struct data*)data)->x;
  double* y = ((struct data*)data)->y;
  fitmode fit_mode =((struct data*)data)->fit_mode;

  for(size_t i=0; i < n; ++i){
    double Yi = fitting_function(x[i], params, fit_mode);
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
void fit_points_with_function(const nvector& interval, const nvector& points, fitting_parameters& fit_parameters, fitmode fit_mode, const unsigned int& verbose){

  const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
  gsl_multifit_nlinear_workspace *w;
  gsl_multifit_nlinear_fdf fdf;
  gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();

  const size_t n = interval.size();
  const size_t p = NUMBER_OF_FITTING_PARAMETERS;

  /* cout of the points to fit
  std::cout << "[";
  for(size_t i = 0; i < interval.size()-1; ++i){
    std::cout << interval[i] << ",";
  }
  std::cout << interval[interval.size()-1] << "]" << std::endl;

  std::cout << "[";
  for(size_t i = 0; i < points.size()-1; ++i){
    std::cout << points[i] << ",";
  }
  std::cout << points[points.size()-1] << "]" << std::endl;
  */

  /* define the function to be minimized */
  gsl_vector* f;
  gsl_matrix* J;
  /* define the covariance matrix */
  gsl_matrix* covar = gsl_matrix_alloc(p,p);

  double x[n], y[n], weights[n], a_init[p];
  nvector initial_guess = guess_initial_fit_parameters(interval, points, fit_mode);

  /* this is the data to be fitted */
  for (size_t i=0; i < n; ++i){
    x[i] = double(interval[i]);
    y[i] = double(points[i]);
    weights[i] = 1.;
  }

  /* set the initial guess for the fitting parameters */
  for(size_t i = 0; i < p; ++i){
    a_init[i] = double(initial_guess[i]);
  }

  if(verbose > 0){
    std::cout << "Starting the fit with the following estimates on the parameters : ";
    for(size_t i=0; i < p; ++i){
      std::cout << a_init[i] << " ";
    }
    std::cout << std::endl;
  }

  struct data d = {p, n, x, y, fit_mode};

  gsl_vector_view a = gsl_vector_view_array(a_init, p);
  gsl_vector_view wts = gsl_vector_view_array(weights, n);
  gsl_rng* r;
  double chisq, chisq0;
  int status, info;

  const double atol = 1e-8;
  const double gtol = 1e-8;
  const double ftol = 0.0;

  gsl_rng_env_setup();
  r = gsl_rng_alloc(gsl_rng_default);


  /* define the function to be minimized*/
  fdf.f = function_to_fit;
  fdf.df = function_to_fit_df;
  fdf.fvv = NULL;
  fdf.n = n;
  fdf.p = p;
  fdf.params = &d;


  /* allocate workspace with default parameters */
  w = gsl_multifit_nlinear_alloc(T, &fdf_params, n, p);

  /* initialize solver with starting points and weights */
  gsl_multifit_nlinear_winit(&a.vector, &wts.vector, &fdf, w);

  /* compute initial cost function
  f = gsl_multifit_nlinear_residual(w);
  gsl_blas_ddot(f,f, &chisq0);
  */

  /* solve the system with a maximum of 1000 iterations */
  status = gsl_multifit_nlinear_driver(1000, atol, gtol, ftol, NULL, NULL, &info, w);

  /* compute covariance of best fit parameters */
  J = gsl_multifit_nlinear_jac(w);
  gsl_multifit_nlinear_covar(J, 0.0, covar);

  /* compute final cost
  double chisq;
  gsl_blas_ddot(f,f,&chisq);
  */

  for(size_t i=0; i < w->x->size; ++i){
    gsl_vector_set(fit_parameters.fit_parameters, i, gsl_vector_get(w->x,i));
    gsl_vector_set(fit_parameters.error, i, sqrt(gsl_matrix_get(covar,i,i)));
  }

  gsl_multifit_nlinear_free(w);
  gsl_rng_free(r);

  return;
}

nvector guess_initial_fit_parameters(const nvector& x, const nvector& y, fitmode fit_mode){
  const size_t p = NUMBER_OF_FITTING_PARAMETERS;
  nvector a_init;

  switch(fit_mode){
    case sigmoidal: {
      /* for the estimate on the solution, we simply take the x whose y is closer to 0 */
      ntype approx_sol = x[0], val_of_approx_sol = y[0];
      for(size_t i = 1; i < x.size(); ++i){
        /* if y[i] is closer to zero than the previous value, change the estimate */
        if(y[i]*y[i] < val_of_approx_sol*val_of_approx_sol){
          approx_sol = x[i];
          val_of_approx_sol = y[i];
        }
      }

      a_init.push_back(approx_sol);

      /* for the estimate on the scale, we take 4/(interval length), (it's an okay approx)*/
      ntype min_el = x[0], max_el = x[x.size()-1];
      for(size_t i=0; i < x.size(); ++i){
        if(x[i] < min_el){
          min_el = x[i];
        }
        if(x[i] > max_el){
          max_el = x[i];
        }
      }
      a_init.push_back(4/(max_el-min_el));
      break;
    }
    default:{
      std::cerr << "This type of fitting has not been implemented yet or does not exist"<<std::endl;
      std::cerr << "Aborting the simulation now"<<std::endl;
      abort();
      break;
    }
  }

  return a_init;
}

statistics estimate_delta_crit_from_interval(const nvector& interval, const nvector& extinctions, const Metaparameters& m, delta_solver delta_solver){
  double x_lo = interval[0];
  double x_hi = interval[interval.size()-1];

  eqmode equilibrium = delta_solver.eq_mode;

  /*so we have this interval of x-points with their y-values, we choose to make a
   polynomial fit around it and take delta critical as the root of that polynomial */

  /* first select the points to fit */
  nvector x_points_to_fit = interval;
  nvector y_points_to_fit = extinctions;

  std::string fit_type;

  switch(delta_solver.fit_mode){
    case sigmoidal:
      NUMBER_OF_FITTING_PARAMETERS = 2;
      fit_type = "sigmoidal";
      break;
    case polynomial:
      NUMBER_OF_FITTING_PARAMETERS = 4;
      fit_type = "degree 3 polynomial";
      break;
    default:
      std::cerr << "This type of fitting has not been implemented yet or does not exist." << std::endl;
      std::cerr << "The program will now abort" << std::endl;
      abort();
      break;
  }

  /* Then we actually find the parameters that fit our choice of function best */
  unsigned int number_of_fitting_parameters = NUMBER_OF_FITTING_PARAMETERS;

  /* define fitting parameters */
  gsl_vector* fit_parameters = gsl_vector_alloc(number_of_fitting_parameters);
  gsl_vector* error = gsl_vector_alloc(number_of_fitting_parameters);

  fitting_parameters fitting_parameters = {fit_parameters, error};
  if(m.verbose > 0){
    std::cout << "Now fitting the " << x_points_to_fit.size() << " points chosen into the specific function (";
    std::cout << fit_type << ")"<< std::endl;
  }
  fit_points_with_function(x_points_to_fit, y_points_to_fit, fitting_parameters, delta_solver.fit_mode, m.verbose);

  /* with the fitting parameters estimated, we can actually solve for Delta numerically */
  statistics delta_crit = solve_for_delta_with_fit(fitting_parameters, x_lo, x_hi, m, delta_solver);
  gsl_vector_free(fit_parameters);
  if(m.verbose > 0){
    std::cout << "Zero estimated at " << delta_crit.mean << "+/-" << delta_crit.std_deviation  << std::endl;
  }

  /* finally, we return the estimated value */
  return delta_crit;
}
