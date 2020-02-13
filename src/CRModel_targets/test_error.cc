#include "CRModel.h"

/*  This script quantifies how many feasible systems you build for certain parameters. This will allow us
    to increase syntrophy while still building a reasonable amount of systems for the same time */
int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    CRModel model1(metaparams);
    std::cout << "Model 1 : " << &model1 << std::endl;
    std::cout << "  Model parameters : " << model1.get_model_parameters() << std::endl;
    std::cout << "    Parameter set : " << model1.get_parameter_set() << std::endl;
    CRModel model2(model1);
    std::cout << "Model 2: " << &model2 << std::endl;
    std::cout << "  Model parameters : " << model2.get_model_parameters() << std::endl;
    std::cout << "    Parameter set : " << model2.get_parameter_set() << std::endl;


    // CRModel model3(model2);
    // std::cout << "Model 3: " << &model3 << std::endl;
    // std::cout << "  Model parameters : " << model3.get_model_parameters() << std::endl;
    // std::cout << "  Equilibria : " << model3.get_equilibrium_abundances() << std::endl;


  }catch(error e){
    e.handle();
  }
  return 0;
}
