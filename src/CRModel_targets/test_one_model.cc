#include <iostream>
#include <fstream>
#include "CRModel.h"

using namespace std;

/*  This script quantifies how many feasible systems you build for certain parameters. */
int main(int argc, char * argv[]){
  try{

    // load metaparameters
    Metaparameters metaparams(argc, argv);

    // create model
    CRModel test_model(metaparams);

    // display alpha matrix
    nmatrix alpha = test_model.get_model_parameters()->get_parameter_set().alpha;
    std::cout << alpha << std::endl;


  }catch(error e){
    e.handle();
  }
  return 0;
}
