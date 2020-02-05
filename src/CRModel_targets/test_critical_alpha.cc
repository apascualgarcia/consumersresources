#include <iostream>
#include <fstream>
#include "CRModel.h"

using namespace std;

/*  This script quantifies how many feasible systems you build for certain parameters. This will allow us
    to increase syntrophy while still building a reasonable amount of systems for the same time */
int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::cout << "Critical alpha for this matrix : " << compute_critical_alpha(metaparams, 1e-4, fitmode(sigmoidal)) << std::endl;
  }catch(error e){
    e.handle();
  }

  return 0;
}
