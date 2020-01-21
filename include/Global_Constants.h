#ifndef GLOBAL_CONSTANTS_H
#define GLOBAL_CONSTANTS_H


const unsigned int print_precision = 7;
const ntype EIGENSOLVER_PRECISION = 1e-15;

/*  At each step the solver changes the step size such that the error level for each component is
    D_i = INTEGRATOR_ABS_PRECISION + INTEGRATOR_RELATIVE_PRECISION *(a_y |y_i| + adydt h |y'i|) */
//const double INTEGRATOR_ABS_PRECISION = std::numeric_limits<double>::epsilon();
const double INTEGRATOR_ABS_PRECISION = 1e-6;
const double INTEGRATOR_REL_PRECISION = 0;
const int INDICES_FOR_AVERAGE=10;

extern std::mt19937 random_engine;

#endif
