#include "../../include/CRModel.h"

int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::cout << find_feasability_probability(metaparams) << std::endl;

  }catch(error e){
    e.handle();
  }

  return 0;
}
