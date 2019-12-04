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
#include "Functions/structural_stability.h"
#include "Functions/additional_functions.h"
#include "Functions/input_reading.h"
#include "Functions/model_characteristics.h"
#include "Functions/equilibrium_study.h"


const unsigned int print_precision = 4;
const ntype EIGENSOLVER_PRECISION = 1e-15;

extern std::mt19937 random_engine;

#endif
