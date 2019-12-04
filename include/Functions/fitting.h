#ifndef FITTING_H
#define FITTING_H

#include <gsl/gsl_vector.h>

#define NUMBER_OF_FITTING_PARAMETERS 4;
double choice_of_fitting_function(double x, void* params);
int function_to_fit(const gsl_vector* , void* , gsl_vector*);


#endif
