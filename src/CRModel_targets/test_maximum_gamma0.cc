#include "../../include/CRModel.h"

int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    nmatrix points = metaparams.load_volume();
    std::cout << "Points:" << points << std::endl;

  }catch(error e){
    e.handle();
  }

  return 0;
}
