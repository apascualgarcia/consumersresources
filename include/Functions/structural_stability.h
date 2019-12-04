#ifndef STRUCTURAL_STABILITY_H
#define STRUCTURAL_STABILITY_H

#include "../Classes/Custom_types.h"
#include "../Classes/Metaparameters.h"

double compute_critical_Delta(Metaparameters, ntype);
double compute_critical_Delta(Metaparameters, ntype, eqmode);

double estimate_delta_crit_from_interval(const nvector&, const nvector&, const Metaparameters&,eqmode);

void write_av_number_extinctions_delta_interval(Metaparameters* , const nvector& , unsigned int Nsimul = 500);
void write_prob_greater_than_one_delta_interval(Metaparameters*, const nvector&, unsigned int Nsimul=500);

#endif
