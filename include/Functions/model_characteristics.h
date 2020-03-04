#ifndef MODEL_CHARACTERISTICS_H
#define MODEL_CHARACTERISTICS_H


int ode_equations_of_evolution(double, const double[], double[], void*);
/* those are the effective equations of evolution, i.e. with resources assumed slowly evolving*/
int effective_ode_equations_of_evolution(double t, const double y[], double f[], void* params);

#endif
