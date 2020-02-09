#ifndef GLOBAL_CONSTANTS_H
#define GLOBAL_CONSTANTS_H

const unsigned int print_precision = 15;
const ntype EIGENSOLVER_PRECISION = 1e-15;

/*  At each step the solver changes the step size such that the error level for each component is
    D_i = INTEGRATOR_ABS_PRECISION + INTEGRATOR_RELATIVE_PRECISION *(a_y |y_i| + adydt h |y'i|) */
//const double INTEGRATOR_ABS_PRECISION = std::numeric_limits<double>::epsilon();
const double INTEGRATOR_ABS_PRECISION = 1e-20;
const double INTEGRATOR_REL_PRECISION = 1e-10;
const double INTEGRATOR_ZERO=0.;
const int INDICES_FOR_AVERAGE=10;

/*  in our unit system, we say a species is extinct if it's abundance goes below
    the value of this constant */
const double SPECIES_EXTINCT = 1e-8;

const unsigned int COMPUTE_CRITICAL_DELTA_POLY_FIRST_RUN=50;
const unsigned int COMPUTE_CRITICAL_DELTA_POLY_SECOND_RUN=100;
const unsigned int COMPUTE_CRITICAL_DELTA_POLY_INTERVAL_LENGTH=10;

const unsigned int COMPUTE_CRITICAL_DELTA_SIG_FIRST_RUN=50;
const unsigned int COMPUTE_CRITICAL_DELTA_SIG_SECOND_RUN=200;
const unsigned int COMPUTE_CRITICAL_DELTA_SIG_INTERVAL_LENGTH=50;

const ntype NUMERICAL_ERROR = std::nanl("1");

extern std::mt19937 random_engine;

#endif
