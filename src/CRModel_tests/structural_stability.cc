#include <iostream>
#include "CRModel.h"
using namespace std;

int main(int argc, char * argv[]){
  Metaparameters metaparams(argc, argv);
  initialize_random_engine(metaparams);
  compute_average_extinction(&metaparams, 0.032, 500);

  return 0;
}
