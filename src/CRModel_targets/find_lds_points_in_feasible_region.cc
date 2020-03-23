#include <iostream>
#include <fstream>
#include "CRModel.h"

using namespace std;

/*  This script quantifies how many feasible systems you build for certain parameters. */
int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::vector<std::string> matrix_list = load_food_matrix_list(metaparams.foodmatrixpath);
    std::ofstream myfile = open_external_file_truncate(metaparams.save_path);
    std::string syntrophy_folder=metaparams.syntrophy_matrix_path;
    nmatrix common_feasible_volume=metaparams.load_volume();
    unsigned int Nsimuls=100;
    if(metaparams.verbose>0){
      std::cout<< common_feasible_volume.size() << " points in cfv " << metaparams.volume_of_interest_path << std::endl;
    }
    for(auto mat: matrix_list){
      metaparams.foodmatrixpath=mat;
      myfile << metaparams.foodmatrixpath << " ";
      metaparams.syntrophy_matrix_path=optimal_alpha_matrix_path_from_syntrophy_folder(metaparams);
      for(auto point : common_feasible_volume){
        metaparams.gamma0=point[0];
        metaparams.S0=point[1];
        std::cout << point << std::endl;
        ntype proba_lds=0.;
        for(size_t i=0; i < Nsimuls; ++i){
          CRModel model(metaparams);
          if(model.is_dynamically_stable()){
            proba_lds+=(1./Nsimuls);
          }
        }
        myfile << metaparams.gamma0 << " " << metaparams.S0 << " " << proba_lds << " ";
      }
      myfile << std::endl;
      if(metaparams.verbose>0){
        std::cout << "Wrote lds points for matrix " << metaparams.foodmatrixpath << std::endl;
      }
      metaparams.syntrophy_matrix_path=syntrophy_folder;
    }
    myfile.close();

  }catch(error e){
    e.handle();
  }
  return 0;
}
