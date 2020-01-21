#include "CRModel.h"

/*  This script quantifies how many feasible systems you build for certain parameters. This will allow us
    to increase syntrophy while still building a reasonable amount of systems for the same time */
int main(int argc, char * argv[]){
  Metaparameters metaparams(argc, argv);

  std::ofstream myfile;
  myfile.open(metaparams.save_path, std::ofstream::out|std::ofstream::trunc);
  if(not(myfile.is_open())){
    std::cerr << "Could not open file " << std::endl;
  }else{
    CRModel model(metaparams);
    eqmode equilibre(convergence);
    writemode ecriture(true, myfile);

    model.perturb_parameters(0.05);
    model.evolve_until_equilibrium(INTEGRATOR_ABS_PRECISION, equilibre,ecriture);
    if(metaparams.verbose > 0){
      std::cout << "Finished writing the time evolution until equilibrium" << std::endl;
    }
  }
  myfile.close();


  return 0;
}
