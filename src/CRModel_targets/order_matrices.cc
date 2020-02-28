#include "../../include/CRModel.h"

int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    CRModel model(metaparams);

    nmatrix Gamma = model.get_Gamma_matrix();
    nmatrix Beta = model.get_Beta_matrix();

    std::cout << "Gamma*Beta=" << std::endl << Gamma*Beta << std::endl;

  }catch(error e){
    e.handle();
  }

  return 0;
}
