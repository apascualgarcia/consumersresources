#include "CRModel.h"

/*  This script quantifies how many feasible systems you build for certain parameters. This will allow us
    to increase syntrophy while still building a reasonable amount of systems for the same time */
int main(int argc, char * argv[]){
  Metaparameters metaparams(argc, argv);

  std::ofstream myfile, eff_file;
  myfile.open(metaparams.save_path, std::ofstream::out|std::ofstream::trunc);
  eff_file.open(metaparams.save_path + "_eff", std::ofstream::out|std::ofstream::trunc);
  if(not(myfile.is_open()) or not(eff_file.is_open())){
    std::cerr << "Could not open file " << std::endl;
  }else{
    CRModel model(metaparams);
    eqmode equilibre(convergence);
    writemode ecriture(true, myfile);

    model.perturb_parameters(0.2);
    model.evolve_until_equilibrium(metaparams.convergence_threshold, equilibre,ecriture);
    //std::cout << "Actual eigenvalues at equilibrium : " << model.eigenvalues_at_equilibrium() << std::endl;
    //std::cout << "Effective eigenvalues at equilibrium : " <<  eff_model.eigenvalues_at_equilibrium() << std::endl;

    if(metaparams.verbose > 0){
      std::cout << "Finished writing the time evolution until equilibrium" << std::endl;
    }
  }
  myfile.close();
  eff_file.close();


  return 0;
}
