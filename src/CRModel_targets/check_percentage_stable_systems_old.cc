#include <iostream>
#include <fstream>
#include "CRModel.h"

using namespace std;

/*  This script quantifies how many feasible systems you build for certain parameters. */
int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::vector<std::string> matrix_list = load_food_matrix_list(metaparams.foodmatrixpath);
    nvector gamma0_interval=linear_interval(0.01, 1., 5);
    unsigned int Nsimuls=100;

    std::ofstream myfile = open_external_file_truncate(metaparams.save_path);
    for(auto mat: matrix_list){
      metaparams.foodmatrixpath=mat;
      myfile << metaparams.foodmatrixpath << " ";
      nmatrix common_feasible_volume = metaparams.common_feasible_volume(10);
      for(auto points : common_feasible_volume){
        metaparams.gamma0=points[0];
        metaparams.S0=points[1];
        ntype dynamically_stab=0.;
        for(size_t i=0; i < Nsimuls; ++i){
          CRModel model(metaparams);
          if(model.is_dynamically_stable()){
            dynamically_stab+=1./Nsimuls;
          }
        }
        myfile << points[0] << " " << points[1] << " " << dynamically_stab << " ";
      }
      myfile << std::endl;
    }
    myfile.close();

  }catch(error e){
    e.handle();
  }
  return 0;
}
