#include "../../include/CRModel.h"
#include <iostream>
int main(int argc, char * argv[]){
  Metaparameters metaparams(argc, argv);
  initialize_random_engine(metaparams);

  stability_metrics stab_metrics = compute_stability_metrics(metaparams);
  std::cout << "Stability_metrics : " << std::endl;
  std::cout << " Resilience : " << stab_metrics.resilience.mean << "+/-" << stab_metrics.resilience.std_deviation << std::endl;

  return 0;
}
