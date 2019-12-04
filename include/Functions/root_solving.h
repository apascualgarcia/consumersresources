#ifndef ROOT_SOLVING_H
#define ROOT_SOLVING_H

#include "../Classes/Custom_types.h"
#include "../Classes/Metaparameters.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>


/* finds a rough interval of Npoints around the zero of function f for the parameters s */
nvector find_rough_interval(gsl_function* f, unsigned int Npoints, unsigned int);

double function_av_extinct_solver(double, void*);
double solve_for_delta_with_fit(const gsl_vector*, double&, double&, const Metaparameters&, eqmode);

#endif
