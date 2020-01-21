#ifndef STRUCTURAL_STABILITY_H
#define STRUCTURAL_STABILITY_H

#include "../Classes/Custom_types.h"
#include "../Classes/Metaparameters.h"


statistics compute_critical_Delta(Metaparameters, stabilitymode stab_mode = structural);
statistics compute_critical_Delta(Metaparameters, delta_solver);

statistics compute_critical_alpha(Metaparameters&, ntype, fitmode);

statistics estimate_delta_crit_from_interval(const nvector&, const nvector&, const Metaparameters&,delta_solver);
statistics estimate_alpha_crit_from_interval(const nvector&, const nvector&, const Metaparameters&,fitmode);

/* compute the different metrics for the stability when we perturb the system by delta */
stability_metrics compute_stability_metrics(Metaparameters&, const ntype& delta, unsigned int Nsimul=100, stabilitymode stab_mode=dynamical);
/* computes the proportion of stable, marginally stable and unstable systems */
stability compute_proportion_stability(Metaparameters&, unsigned int Nsimuls);

void write_av_number_extinctions_delta_interval(Metaparameters* , const nvector& , unsigned int Nsimul = 500);
void write_prob_greater_than_one_delta_interval(Metaparameters*, const nvector&, unsigned int Nsimul=500);



#endif
