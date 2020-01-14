#ifndef CRMODEL_H
#define CRMODEL_H

#include "Classes/Custom_types.h"
#include "Classes/Extinction.h"
#include "Classes/ParameterSet.h"
#include "Classes/Solver.h"
#include "Classes/Metaparameters.h"
#include "Classes/Dynamical_variables.h"
#include "Classes/Model_Parameters.h"
#include "Classes/Consumer_Resource_Model.h"

#include "Functions/output_writing.h"
#include "Functions/fitting.h"
#include "Functions/root_solving.h"
#include "Functions/matrices_building.h"
#include "Functions/study_stability.h"
#include "Functions/additional_functions.h"
#include "Functions/input_reading.h"
#include "Functions/model_characteristics.h"
#include "Functions/equilibrium_study.h"


#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>


#include<limits>
const unsigned int print_precision = 4;
const ntype EIGENSOLVER_PRECISION = 1e-15;

/*  At each step the solver changes the step size such that the error level for each component is
    D_i = INTEGRATOR_ABS_PRECISION + INTEGRATOR_RELATIVE_PRECISION *(a_y |y_i| + adydt h |y'i|) */
//const double INTEGRATOR_ABS_PRECISION = std::numeric_limits<double>::epsilon();
const double INTEGRATOR_ABS_PRECISION = 1e-6;
const double INTEGRATOR_REL_PRECISION = 0;

extern std::mt19937 random_engine;

#endif
