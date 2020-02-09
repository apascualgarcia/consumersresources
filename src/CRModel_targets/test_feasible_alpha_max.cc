#include "../../include/CRModel.h"

int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    ntype eps_min=0.1, eps_max=0.25;
    unsigned int eps_points=10;

    nvector eps_interval = linear_interval(eps_min, eps_max, eps_points);
    for(size_t i=0; i < eps_points; ++i){
      metaparams.epsilon=eps_interval[i];
      std::cout << "For epsilon=" << metaparams.epsilon << ", maximum feasible alpha0 is ";
      std::cout << metaparams.feasible_alpha_max() << std::endl;

    }



  }catch(error e){
    e.handle();
  }
  return 0;
}
