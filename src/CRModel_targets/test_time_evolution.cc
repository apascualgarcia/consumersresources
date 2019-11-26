#include "CRModel.h"
#include <iostream>

int main(int argc, char * argv[]){
  Metaparameters metaparams(argc, argv);
  initialize_random_engine(metaparams);
  for(size_t i = 0 ; i < 1000; ++i){
    CRModel model(metaparams);
    model.perturb_parameters(metaparams.perturb_parameters);
    Extinction extinction = model.evolve_until_equilibrium(1e-6, eqmode(oneextinct));
  }

  return 0;
}
