#include "../../include/CRModel.h"

int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    initialize_random_engine(metaparams);

    CRModel model(metaparams);

    std::cout << "Initial equilibrium : " << std::endl <<model.get_first_equilibrium() << std::endl;
    Extinction extinct = model.evolve_until_equilibrium_from_abundances(model.perturb_abundances(0.5));
    std::cout << "Equilibrium after perturbation : " << std::endl;
    std::cout << extinct.new_Req << std::endl;
    std::cout << extinct.new_Seq << std::endl;
    std::cout << "Time to get convergence : " << extinct.t_eq << std::endl;
  }catch(error e){
    e.handle();
  }
  return 0;
}
