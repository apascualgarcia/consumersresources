#include <iostream>
#include <fstream>
#include "CRModel.h"

using namespace std;

/*  This script quantifies how many feasible systems you build for certain parameters. */
int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    nvector alpha0_interval=linear_interval(0.04, 0.047, 50);
    for(size_t i =0; i < alpha0_interval.size(); ++i){
      metaparams.alpha0 = alpha0_interval[i];
      std::cout << "For alpha0 = " << alpha0_interval[i] << ", feasibility proba: " << find_feasability_probability(metaparams, 1e3) << std::endl;
    }

  }catch(error e){
    e.handle();
  }
  return 0;
}
