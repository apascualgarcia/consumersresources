#ifndef ROOT_SOLVING_H
#define ROOT_SOLVING_H

#include "../Classes/Custom_types.h"
#include "../Classes/Metaparameters.h"

struct delta_solver;

/* finds a rough interval of Npoints around the zero of function f for the parameters s */
nvector find_rough_interval(gsl_function* f, unsigned int Npoints, unsigned int, fitmode, interval);

nvector find_rough_interval_polynomial_fit(gsl_function* f, unsigned int Npoints, unsigned int, interval first_guess);
nvector find_rough_interval_sigmoidal_fit(gsl_function* f, unsigned int Npoints, unsigned int, interval first_guess);

/* will give a point where the function f is roughly 0 (can shift this with the shift parameter for f) */
double find_zero(gsl_function* f, unsigned int verbose, interval);

double function_av_extinct_solver(double, void*);
statistics solve_for_delta_with_fit(fitting_parameters&, double&, double&, const Metaparameters&, delta_solver);

/* takes an alpha0 and Solver Parameters as input, gives back feasability */
double function_proba_feasability_solver(double, void*);
/* takes an alpha0 and Solver Parameters as input, gives back change of locally dynamically stable */
double function_proba_dynamical_stability_solver(double alpha0, void*);
/* takes an S0 and Solver Parameters as input, gives back feasability */
double function_proba_feasability_solver_S0(double S, void* params);
/* takes a gamma0 and Solver Parameters as input, gives back feasability */
double function_proba_feasability_solver_gamma0(double gamma, void* params);



#endif
