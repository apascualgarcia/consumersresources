#include <iostream>
#include <fstream>
#include "CRModel.h"

using namespace std;

/*  This script quantifies how many feasible systems you build for certain parameters. */
int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::vector<std::string> matrix_list = load_food_matrix_list(metaparams.foodmatrixpath);
    std::string syntrophy_folder=metaparams.syntrophy_matrix_path;
    std::vector<alphamode> alpha_modes = {random_structure, no_release_when_eat, optimal_matrix};
    ntype max_alpha0=0.;

    nvector g0_interval=linear_interval(0.01, 1, 30);
    nvector S0_interval=linear_interval(0.01, 1, 30);

    for(auto g : g0_interval){
      for(auto S : S0_interval){
        metaparams.gamma0=g;
        metaparams.S0=S;
        for(auto alpha: alpha_modes){
          metaparams.alpha_mode=alpha;
          for(auto mat : matrix_list){
            metaparams.foodmatrixpath=mat;
            metaparams.syntrophy_matrix_path=optimal_alpha_matrix_path_from_syntrophy_folder(metaparams);
            ntype local_alpha_max = metaparams.feasible_alpha_max(1e-5);
            if(local_alpha_max > max_alpha0){
              max_alpha0 = local_alpha_max;
            }
            metaparams.syntrophy_matrix_path=syntrophy_folder;
          }
        }
      }
    }

    std::cout << "Maximum feasible alpha0 (all modes) = " << max_alpha0 << std::endl;

  }catch(error e){
    e.handle();
  }
  return 0;
}
