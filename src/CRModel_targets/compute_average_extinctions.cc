#include <iostream>
#include <fstream>
#include "CRModel.h"

using namespace std;

int main(int argc, char * argv[]){
  Metaparameters metaparams(argc, argv);
  initialize_random_engine(metaparams);
  std::ofstream myfile;
  std::string spath = metaparams.save_path;
  myfile.open(spath, std::ofstream::out | std::ofstream::trunc);
  if(not(myfile.is_open())){
    std::cerr << "Could not open " << spath << " for writing the extinctions " << std::endl;
  }else{
    int Nsimul = 10000;
    for(size_t i = 0 ; i < Nsimul; ++i){
      CRModel model(metaparams);
      model.perturb_parameters(metaparams.perturb_parameters);
      Extinction ext = model.evolve_until_equilibrium(1e-6);
      myfile  << ext.extinct << " " << ext.t_eq << std::endl;
    }
    myfile << std::endl;
  }
  myfile.close();
  return 0;
}
