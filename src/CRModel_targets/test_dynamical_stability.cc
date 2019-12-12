#include "../../include/CRModel.h"

int main(int argc, char * argv[]){
  Metaparameters metaparams(argc, argv);
  initialize_random_engine(metaparams);

  CRModel model(metaparams);
  std::cout << "Initial equilibrium : " << std::endl <<model.get_first_equilibrium() << std::endl;
  Extinction extinct = model.evolve_until_equilibrium_from_abundances(model.perturb_abundances(0.));
  std::cout << "Equilibrium after perturbation : " << std::endl;
  std::cout << extinct.new_Req << std::endl;
  std::cout << extinct.new_Seq << std::endl;
  std::cout << "Distance between two equilibria : " << distance_between_equilibria(extinct)<<std::endl;

  return 0;
}
