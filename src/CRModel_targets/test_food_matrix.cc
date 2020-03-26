#include "../../include/CRModel.h"

int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    metaparams.NS=25;
    metaparams.NR=50;
    nmatrix gamma = load_food_matrix(metaparams);
    metaparams.NS=50;
    metaparams.NR=25;
    metaparams.foodmatrixpath = "data_output/error_matrix.out";
    nmatrix alpha = load_food_matrix(metaparams);
    display_food_matrix(std::cout, alpha);
    std::cout << std::endl;
    for(size_t i=0; i < 1000; ++i){
      modify_column(alpha, gamma, true);
    }




  }catch(error e){
    e.handle();
  }

  return 0;
}
