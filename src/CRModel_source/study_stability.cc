#include "../../include/CRModel.h"
#include <gsl/gsl_statistics.h>

stability_metrics compute_stability_metrics(Metaparameters& m, unsigned int Nsimul){
  nvector resiliences;
  for(size_t i = 0; i < Nsimul; ++i){
    CRModel model(m);
    resiliences.push_back(model.get_resilience_dynamical_stability());
  }
  ntype av_resilience = mean(resiliences);
  ntype std_resilience = standard_dev(resiliences);
  /* /!\ CAREFUL resiliences is shuffled in the gsl_stats_median part */
  ntype median_resilience = 0;

  statistics resilience = {av_resilience, std_resilience, median_resilience};
  stability_metrics stab_metr = {resilience};
  return stab_metr;
}
