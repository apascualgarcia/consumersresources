#include <iostream>
#include "CRModel.h"
#include <fstream>
using namespace std;

int main(int argc, char * argv[]){
    Metaparameters metaparams(argc, argv);
    initialize_random_engine(metaparams);
    std::vector<std::string> matrices_path=load_food_matrix_list(metaparams.foodmatrixpath);
    std::string syntrophy_folder=metaparams.syntrophy_matrix_path;

    /* CAREFUL WE SOLVE FOR DYNAMICAL STABILITY HERE */
    stabilitymode stab = dynamical;

    std::ofstream myfile;
    myfile.open(metaparams.save_path, std::ios::app);
    if(metaparams.verbose > 0){
      std::cout << "The critical deltas for the following matrices will be computed :";
      for(size_t i = 0 ; i < matrices_path.size(); ++i){
        std::cout << std::endl << matrices_path[i];
      }
      std::cout << std::endl;
    }
    bool save_success(false);
    if(not(myfile.is_open())){
      std::cerr << "Could not open " << metaparams.save_path << " to write the new equilibrium of the system" << std::endl;
    }else{
      if(metaparams.verbose > 0){
        std::cout << "Successfully opened " << metaparams.save_path <<", attempting now to find the critical delta of every listed matrix " << std::endl;
      }
      save_success = true;
      for(size_t i = 0; i < matrices_path.size();++i){
          metaparams.foodmatrixpath = matrices_path[i];
          metaparams.syntrophy_matrix_path=optimal_alpha_matrix_path_from_syntrophy_folder(metaparams);
          std::cout << "Feasability probability is " << find_feasability_probability(metaparams) << std::endl;
          statistics delta = compute_critical_Delta(metaparams, stab);
          std::cout << "Computed critical delta for " << matrices_path[i] << std::endl;
          myfile << matrices_path[i] << " " << delta.mean_ << " " << delta.std_deviation_ << std::endl;
          metaparams.syntrophy_matrix_path=syntrophy_folder;
      }
    }
    myfile.close();
  return 0;
}
