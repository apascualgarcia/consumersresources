#include <iostream>
#include "CRModel.h"

using namespace std;

int main(int argc, char * argv[]){
  Metaparameters metaparams(argc, argv);
  initialize_random_engine(metaparams);
  foodmatrix food_matrix = load_food_matrix(metaparams);
  CRModel model(food_matrix, metaparams);
  cout << model << endl;

  return 0;
}
