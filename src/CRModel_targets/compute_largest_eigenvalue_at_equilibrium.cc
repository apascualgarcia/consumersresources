#include "CRModel.h"

int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::vector<std::string> matrices=load_food_matrix_list(metaparams.foodmatrixpath);
    std::string syntrophy_folder=metaparams.syntrophy_matrix_path;

    std::ofstream myfile = open_external_file_truncate(metaparams.save_path);
    unsigned int Npoints(100), Nsimuls(1000);

    /*  Algorithmic procedure : for each matrix we find the largest feasible alpha0
        we then compute many points to see where lies the transition from dynamically
        stable to unstable */
    for(size_t i=0; i < matrices.size(); ++i){
      metaparams.foodmatrixpath = matrices[i];
      metaparams.syntrophy_matrix_path=optimal_alpha_matrix_path_from_syntrophy_folder(metaparams);
      if(metaparams.verbose>0){
        std::cout << "Computing points for consumption matrix " << metaparams.foodmatrixpath << std::endl;
        std::cout << "\t and syntrophy matrix " << metaparams.syntrophy_matrix_path << std::endl;
      }

      ntype max_alpha0 = metaparams.feasible_alpha_max(1e-6);
      nvector alpha_range = linear_interval(0., max_alpha0, Npoints);
      for(size_t j=0; j < Npoints; ++j){
        metaparams.alpha0 = alpha_range[j];
        CRModel  model(metaparams);
        nctype largest_eigenvalue_at_equilibrium=model.largest_eigenvalue_at_equilibrium();
        for(size_t k=0; k < Nsimuls; ++k){
          CRModel model_test(metaparams);
          nctype local_largest = model_test.largest_eigenvalue_at_equilibrium();
          if(local_largest > largest_eigenvalue_at_equilibrium){
            largest_eigenvalue_at_equilibrium = local_largest;
          }
        }
        myfile << metaparams.alpha0 << " " << largest_eigenvalue_at_equilibrium << " ";
      }
      myfile << std::endl;
      metaparams.syntrophy_matrix_path=syntrophy_folder;
    }
    myfile.close();

  }catch(error e){
    e.handle();
  }

  return 0;
}
