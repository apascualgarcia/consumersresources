#include "../../include/CRModel.h"

int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::vector<std::string> matrix_list=load_food_matrix_list(metaparams.foodmatrixpath);

    nvector expected_S0;
    for(auto mat : matrix_list){
      metaparams.foodmatrixpath=mat;
      foodmatrix gamma = load_food_matrix(metaparams);
      expected_S0.push_back(1./maximum(columns_degrees(gamma)));
      std::cout << "Matrix nest=" << nestedness(gamma) << ", conn=" << connectance(gamma) << " expected coeff=" << 1./maximum(columns_degrees(gamma))<<std::endl;
    }
    std::cout << "Expected S0 = [";
    for(auto S : expected_S0){
      std::cout << S << ",";
    }
    std::cout << "]" << std::endl;
    std::cout << "Expected coefficient : " << metaparams.l0*minimum(expected_S0) << std::endl;

  }catch(error e){
    e.handle();
  }

  return 0;
}
