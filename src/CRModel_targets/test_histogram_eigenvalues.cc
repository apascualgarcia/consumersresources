#include "../../include/CRModel.h"

int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    unsigned int Nsimul=10000;
    std::ofstream file=open_external_file_truncate(metaparams.save_path);
    for(size_t i=0; i < Nsimul; ++i){
      CRModel model(metaparams);
      file << model.largest_eigenvalue_at_equilibrium() <<" ";
    }
    file.close();
  }catch(error e){
    e.handle();
  }

  return 0;
}
