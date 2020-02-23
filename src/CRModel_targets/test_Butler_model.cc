#include "../../include/CRModel.h"


int main(int argc, char * argv[]){
  Metaparameters metaparams(argc, argv);
  ButlerModel model(metaparams);
  std::cout << "Model from Butler's paper: " << std::endl << model << std::endl;
  std::cout << "Largest eigenvalue : " << model.largest_eigenvalue_at_equilibrium() << std::endl;


  return 0;
}
