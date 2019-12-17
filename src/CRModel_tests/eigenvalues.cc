#include <iostream>
#include "CRModel.h"

using namespace std;
std::mt19937 random_device(0);

int main(int argc, char* argv[]){
  Metaparameters metaparams(argc, argv);
  initialize_random_engine(metaparams);
  foodmatrix food_matrix = load_food_matrix(metaparams);
  CRModel model(food_matrix, metaparams);
  model.save_simulation();

  return 0;
}
