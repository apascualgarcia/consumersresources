#include <iostream>
#include "CRModel.h"

using namespace std;

int main(int argc, char * argv[]){
  Metaparameters metaparams(argc, argv);
  initialize_random_engine(metaparams);
  foodmatrix food_matrix = load_food_matrix(metaparams);
  CRModel model(food_matrix, metaparams);
  model.perturb_parameters();
  Extinction extinct = model.evolve_until_equilibrium(1e-6);
  model.save_new_equilibrium(extinct);
  std::cout << "Hello World" << std::endl;

  return 0;
}
