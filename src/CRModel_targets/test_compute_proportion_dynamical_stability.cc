#include "../../include/CRModel.h"

int main(int argc, char * argv[]){
  Metaparameters metaparams(argc, argv);
  unsigned int Nsimuls = 1000;

  stability stab = compute_proportion_stability(metaparams, Nsimuls);
  std::cout << "For the following metaparameters : " << std::endl;
  std::cout << metaparams << std::endl;
  std::cout << "We found :" << std::endl;
  std::cout << "\t Unstable systems : " << 100*stab.unstable << " \%" << std::endl;
  std::cout << "\t Marginally stable systems : " << 100*stab.marginally_stable << " \%" << std::endl;
  std::cout << "\t Stable systems : " << 100*stab.stable << " \%" << std::endl;

  return 0;
}
