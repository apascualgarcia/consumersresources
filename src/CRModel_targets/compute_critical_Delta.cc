#include <iostream>
#include "CRModel.h"

using namespace std;

int main(int argc, char * argv[]){
  Metaparameters metaparams(argc, argv);
  initialize_random_engine(metaparams);
  double critical_delta = compute_critical_Delta(metaparams,0.001);
  std::cout << "Critical delta found for this matrix : " << critical_delta << std::endl;
  return 0;
}
