#include "../../include/CRModel.h"


int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    metaparams.building_mode=use_m;
    metaparams.m0=0.;
    metaparams.alpha0=0.;

    unsigned int Nsimuls = 1000;

    std::ofstream myfile = open_external_file_truncate(metaparams.save_path);
    myfile << metaparams.foodmatrixpath << " ";
    for(size_t i=0; i < Nsimuls; ++i){
      CRModel model(metaparams);
      ncvector eigvals = model.eigenvalues_at_equilibrium();
      myfile << eigvals << " ";
    }

    myfile.close();


  }catch(error e){
    e.handle();
  }
  return 0;
}
