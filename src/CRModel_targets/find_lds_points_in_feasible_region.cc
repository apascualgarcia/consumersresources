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

    for(auto mat: matrix_list){
      metaparams.foodmatrixpath=mat;
      myfile << metaparams.foodmatrixpath << " ";
      metaparams.syntrophy_matrix_path=optimal_alpha_matrix_path_from_syntrophy_folder(metaparams);
      nmatrix lds_points = metaparams.set_of_feasible_points();
      for(auto point : lds_points){
        myfile << point;
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
