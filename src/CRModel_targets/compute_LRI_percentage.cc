#include "../../include/CRModel.h"

int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::vector<std::string> matrices = load_food_matrix_list(metaparams.foodmatrixpath);
    std::string syntrophy_folder = metaparams.syntrophy_matrix_path;

    /* take smallest commun feasible alpha0 */
    metaparams.alpha_mode=optimal_matrix;

    /* simulate Nsimuls systems for each matrix and compute the proportion that is in an LRI regime */
    unsigned int Nsimuls=10000;
    for(size_t i=0; i < matrices.size(); ++i){
      metaparams.foodmatrixpath=matrices[i];
      ntype alpha0max=metaparams.feasible_alpha_max(1e-4);
      metaparams.alpha0=alpha0max;
      std::cout << "Maximum feasible alpha0 " << alpha0max << std::endl;

      //metaparams.syntrophy_matrix_path=optimal_alpha_matrix_path_from_syntrophy_folder(metaparams);

      ntype proba_LRI=0.;
      for(size_t j=0; j < Nsimuls; ++j){
        CRModel model(metaparams);
        if(model.is_in_low_intra_resource_interaction()&&model.is_dynamically_stable()){
          proba_LRI+=1./Nsimuls;
        }
      }
      std::cout << "For this matrix proba of LRI is " << proba_LRI << std::endl;
      metaparams.syntrophy_matrix_path=syntrophy_folder;
    }


  }catch(error e){
    e.handle();
  }
  return 0;
}
