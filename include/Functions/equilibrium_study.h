#ifndef EQUILIBRIUM_STUDY_H
#define EQUILIBRIUM_STUDY_H

#include "../Classes/Custom_types.h"
#include "../Classes/Metaparameters.h"

// for a given set of metaparameters, computes the average extinction of the system
Extinction_statistics compute_average_extinction(Metaparameters*, const ntype &, unsigned int);

double average_number_of_extinctions(double , Metaparameters*, unsigned int);
double probability_of_extinction_greather_than_one(Metaparameters*, const ntype& , unsigned int);

#endif
