#include "CRModel.h"

/*  This script quantifies how many feasible systems you build for certain parameters. This will allow us
    to increase syntrophy while still building a reasonable amount of systems for the same time */
int main(int argc, char * argv[]){
  Metaparameters metaparams(argc, argv);
  CRModel model(metaparams);
  std::cout << "Flux for consumers : " << std::endl;
  for(size_t i=0; i < metaparams.NS; ++i){
    std::cout << model.syntrophy_flux_equilibrium_consumer(i) << " (sy) " << model.consumption_intake_flux_equilibrium_consumer(i) << " (cons)" << std::endl;
  }
  std::cout << std::endl;

  std::cout << "Flux for resources : " << std::endl;
  for(size_t mu=0; mu < metaparams.NR; ++mu){
    std::cout << model.syntrophy_flux_equilibrium_resource(mu) << " (sy) " << model.consumption_flux_equilibrium_resource(mu) << " (cons) " << std::endl;
  }
  std::cout << std::endl;

  return 0;
}
