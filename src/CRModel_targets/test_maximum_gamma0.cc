#include "../../include/CRModel.h"

int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::vector<std::string> matrices=load_food_matrix_list(metaparams.foodmatrixpath);

    for(size_t i=0; i < matrices.size(); ++i){
      metaparams.foodmatrixpath = matrices[i];
      std::cout << "Matrix " << matrices[i] << " : ratio is ";
      ntype amax=  metaparams.feasible_alpha_max(1e-4);
      metaparams.alpha0=amax;
      CRModel model(metaparams);
      nmatrix gamma = model.get_parameter_set()->gamma;
      nmatrix alpha = model.get_parameter_set()->alpha;
      ntype thamax = metaparams.sigma0/(1-metaparams.sigma0)*metaparams.l0/(connectance(alpha)*metaparams.NS*metaparams.S0);
      std::cout << amax/thamax << std::endl;
    }

  }catch(error e){
    e.handle();
  }

  return 0;
}
