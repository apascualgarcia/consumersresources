#include <iostream>
#include "CRModel.h"

using namespace std;

/*  This script quantifies how many feasible systems you build for certain parameters. This will allow us
    to increase syntrophy while still building a reasonable amount of systems for the same time */
int main(int argc, char * argv[]){
  Metaparameters metaparams(argc, argv);
  /* we attempt to build a Nsystems systems for each point*/
  unsigned int Nsystems = 1000;
  /* we check for N_alpha values of alpha between alpha_min and alpha_max*/
  unsigned int N_alpha = 10;
  ntype alpha_min = 0.;
  ntype alpha_max = 4.;

  nvector probability_feasability;

  for(size_t i=0; i<Nsystems;++i){

  }



  return 0;
}
