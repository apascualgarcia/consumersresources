#include "../../include/CRModel.h"

int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);

    std::vector<std::string> matrices = load_food_matrix_list(metaparams.foodmatrixpath);
    ntype S0_min = 1e9;

    for(size_t i=0; i < matrices.size(); ++i){
      metaparams.foodmatrixpath=matrices[i];
      ntype local_S0 = metaparams.feasible_S0_max();
      if(local_S0<S0_min){
        S0_min = local_S0;
      }
    }
    std::cout << "Lowest common S0 is " << S0_min << std::endl;


  }catch(error e){
    e.handle();
  }

  return 0;
}
