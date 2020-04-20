#include "../../include/CRModel.h"

int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::ofstream file=open_external_file_truncate(metaparams.save_path);
    unsigned int Nsimul=100000;
    nmatrix jacobian_at_equilibrium(metaparams.NR+metaparams.NS, nvector(metaparams.NR+metaparams.NS, 0.));

    for(size_t i=0; i < Nsimul; ++i){
      CRModel model(metaparams);
      jacobian_at_equilibrium=jacobian_at_equilibrium+model.jacobian_at_equilibrium()/Nsimul;
    }

    file << jacobian_at_equilibrium << std::endl;
    file.close();
  }catch(error e){
    e.handle();
  }

  return 0;
}
