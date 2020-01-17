#ifndef EQUILIBRIUM_STUDY_H
#define EQUILIBRIUM_STUDY_H

#include "../Classes/Custom_types.h"
#include "../Classes/Metaparameters.h"

/* finds the probability that the metaparameters yield a feasible system (probability estimated with Nruns runs) */
ntype find_feasability_probability(Metaparameters& metaparams, unsigned int Nruns=1000);

// for a given set of metaparameters, computes the average extinction of the system
Extinction_statistics compute_average_extinction(Metaparameters*, const ntype &, unsigned int);

double average_number_of_extinctions(double , Metaparameters*, unsigned int);
double probability_of_extinction_greather_than_one(Metaparameters*, const ntype& delta, unsigned int Nsimul, stabilitymode stab_mode=structural);
/* perturbs a system by delta (structurally or dynamically and computes the av. dist between eq) */
statistics distance_between_equilibria(Metaparameters*, const ntype& delta, unsigned int Nsimul, stabilitymode stab_mode);
statistics angle_between_equilibria(Metaparameters*, const ntype& delta, unsigned int Nsimul, stabilitymode stab_mode);

/* returns 0. if we find one system that goes extinct within the given number of simulations, returns 1. otherwise */
double can_find_one_extinction(Metaparameters*, const ntype & delta, unsigned int Nsimul);
/*  returns 0. if we once find a system where there is no extinction, returns 1. if we don't find any system that has 0 extinctions */
double can_find_zero_extinction(Metaparameters*, const ntype & delta, unsigned int Nsimul);
#endif
