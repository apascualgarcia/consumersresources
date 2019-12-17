#ifndef STRUCTURAL_STABILITY_H
#define STRUCTURAL_STABILITY_H

#include "../Classes/Custom_types.h"
#include "../Classes/Metaparameters.h"


statistics compute_critical_Delta(Metaparameters, ntype, stabilitymode stab_mode = structural);
statistics compute_critical_Delta(Metaparameters, ntype, delta_solver);

statistics estimate_delta_crit_from_interval(const nvector&, const nvector&, const Metaparameters&,delta_solver);

/* compute the different metrics for the stability when we perturb the system by delta */
stability_metrics compute_stability_metrics(Metaparameters&, const ntype& delta, unsigned int Nsimul=100, stabilitymode stab_mode=dynamical);

void write_av_number_extinctions_delta_interval(Metaparameters* , const nvector& , unsigned int Nsimul = 500);
void write_prob_greater_than_one_delta_interval(Metaparameters*, const nvector&, unsigned int Nsimul=500);



#endif
