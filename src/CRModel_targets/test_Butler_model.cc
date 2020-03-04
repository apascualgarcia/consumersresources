#include "../../include/CRModel.h"


int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    metaparams.building_mode=use_m;
    metaparams.alpha0=0.;
    metaparams.m0=0.;

    std::vector<std::string> matrices = load_food_matrix_list(metaparams.foodmatrixpath);

    unsigned int Nsimuls=1000;

    for(size_t i=0; i < matrices.size();++i){
      metaparams.foodmatrixpath=matrices[i];
      CRModel model(metaparams);
      nctype largest_eigenvalue_at_equilibrium = model.largest_eigenvalue_at_equilibrium();
      //std::cout << model << std::endl;
      for(size_t k=0; k < Nsimuls-1; ++k){
        CRModel model_test(metaparams);
        nctype local_largest = model_test.largest_eigenvalue_at_equilibrium();
        if(local_largest > largest_eigenvalue_at_equilibrium){
          largest_eigenvalue_at_equilibrium = local_largest;
        }
      }
      std::cout << "Largest eigenvalue at equilibrium : " << largest_eigenvalue_at_equilibrium << std::endl;
    }
  }catch(error e){
    e.handle();
  }

  return 0;
}
