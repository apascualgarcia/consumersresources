#include "CRModel.h"

/*  This script quantifies how many feasible systems you build for certain parameters. This will allow us
    to increase syntrophy while still building a reasonable amount of systems for the same time */
int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);

    double perturbation = 0.2;
    unsigned int Nsimuls=10;
    double proba_when_struct_perturbed = probability_of_extinction_greather_than_one(&metaparams, perturbation, Nsimuls);
    std::cout << "For convergence threshold="<<metaparams.convergence_threshold;
    std::cout << " we observe a probability of extinction of "  << proba_when_struct_perturbed ;
    std::cout << " when we perturb the model with Delta=" << perturbation;
    std::cout << std::endl;
  }catch(error e){
    e.handle();
  }
  return 0;
}
