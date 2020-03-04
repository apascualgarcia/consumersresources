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
#include "Classes/Effective_Consumer_Resource_Model.h"
#include "Classes/ButlerModel.h"

#include "Functions/output_writing.h"
#include "Functions/fitting.h"
#include "Functions/root_solving.h"
#include "Functions/matrices_building.h"
#include "Functions/study_stability.h"
#include "Functions/additional_functions.h"
#include "Functions/input_reading.h"
#include "Functions/model_characteristics.h"
#include "Functions/equilibrium_study.h"
#include "Functions/optimize_alpha_matrix.h"


#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>

#include<limits>
#include<fstream>
#include<iostream>
#include<iomanip>
#include<cmath>
#include<chrono>
#include<ctime>
#include<complex>
#include<list>
#include<algorithm>
#include<Eigen/Dense>
#include "Global_Constants.h"

#endif
