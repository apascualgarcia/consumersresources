#include <iostream>
#include "CRModel.h"
#include <gsl/gsl_math.h>

using namespace std;

double f(double x){
  double root = 0.1;
  return (x-root)*(x-root)*(x-root);
}


int main(int argc, char * argv[]){
  // Metaparameters metaparams(argc, argv);
  // initialize_random_engine(metaparams);
  // foodmatrix food_matrix = load_food_matrix(metaparams);
  // CRModel model(food_matrix, metaparams);
  // model.perturb_parameters();
  // Extinction extinct = model.evolve_until_equilibrium(1e-6);
  // model.save_new_equilibrium(extinct);
  // std::cout << "Hello World" << std::endl;

  nvector x, y;
  ntype a = 0.0, b=1;
  for(size_t i=0; i < 10; ++i){
    x.push_back(0.+i*(b-a)/9);
    y.push_back(f(x[i]));
  }

  std::cout << " x =" << x << std::endl;
  std::cout << " y =" << y << std::endl;

  estimate_delta_crit_from_interval(x,y);

  return 0;
}
