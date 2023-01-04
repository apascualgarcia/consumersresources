#ifndef FITTING_H
#define FITTING_H

#include "../CRModel.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit_nlinear.h>


extern unsigned int NUMBER_OF_FITTING_PARAMETERS;

struct data {
  size_t p; // number of fitting parameters
  size_t n; // number of points to fit
  double* x;
  double* y;
  fitmode fit_mode;
};

struct fitting_parameters{
  gsl_vector* fit_parameters;
  gsl_vector* error;
};

double polynomial_fit(double, const gsl_vector*);
double sigmoidal_fit(double, const gsl_vector*);
double sigmoidal_fit_erf(double, const gsl_vector*);
double fitting_function(double, const gsl_vector*, fitmode);

nvector guess_initial_fit_parameters(const nvector&, const nvector&, fitmode);

int function_to_fit(const gsl_vector* , void* , gsl_vector*);
int function_to_fit_df(const gsl_vector*, void*, gsl_matrix*);
void callback(const size_t, void*, const gsl_multifit_nlinear_workspace);
void fit_points_with_function(const nvector&, const nvector&, fitting_parameters&, fitmode, const unsigned int& = 1);

void set_jacobian_matrix_fit(const gsl_vector*, void*, gsl_matrix*);
void set_jacobian_matrix_sigmoidal_fit(const gsl_vector*, void*, gsl_matrix*);
void set_jacobian_matrix_polynomial_fit(const gsl_vector*, void*, gsl_matrix*);


#endif
