#include "CRModel.h"

using namespace std;

/*  This script quantifies how many feasible systems you build for certain parameters. */
int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::vector<std::string> matrix_list = load_food_matrix_list(metaparams.foodmatrixpath);
    std::string syntrophy_folder=metaparams.syntrophy_matrix_path;
    std::vector<alphamode> alpha_modes = {optimal_matrix};
    std::ofstream myfile=open_external_file_truncate(metaparams.save_path);
    myfile << "# alpha modes for";
    for(auto alpha : alpha_modes){
      myfile << " " << alpha << ";";
    }
    myfile << std::endl;
    for(auto mat: matrix_list){
      metaparams.foodmatrixpath=mat;
      myfile << metaparams.foodmatrixpath;
      metaparams.syntrophy_matrix_path=optimal_alpha_matrix_path_from_syntrophy_folder(metaparams);
      for(auto alpha : alpha_modes){
        metaparams.alpha_mode=alpha;
        ntype local_alpha_max = metaparams.dynamical_alpha_max(1e-5);
        myfile << " " << local_alpha_max;
        std::cout << "Critical dynamical alpha0 for " << metaparams.foodmatrixpath << " (syntrophy mode "<< alpha << ") is " << local_alpha_max << std::endl;
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
