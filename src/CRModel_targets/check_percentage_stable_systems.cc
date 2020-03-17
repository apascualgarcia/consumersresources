#include <iostream>
#include <fstream>
#include "CRModel.h"

using namespace std;

/*  This script quantifies how many feasible systems you build for certain parameters. */
int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::vector<std::string> matrix_list = load_food_matrix_list(metaparams.foodmatrixpath);
    nvector gamma0_interval=linear_interval(0.01, 1., 100);
    unsigned int Nsimuls=1000;

    std::ofstream myfile = open_external_file_truncate(metaparams.save_path);
    for(auto mat: matrix_list){
      metaparams.foodmatrixpath=mat;
      myfile << metaparams.foodmatrixpath << " ";
      for(auto g : gamma0_interval){
        // S0 only goes up to its maximum value
        nvector S0_interval=linear_interval(0.01, 0.04261718*g-0.00456834, 100);
        for(auto S : S0_interval){
          metaparams.gamma0=g;
          metaparams.S0=S;
          ntype dynamically_stab=0.;
          for(size_t i=0; i < Nsimuls; ++i){
            CRModel model(metaparams);
            if(model.is_dynamically_stable()){
              dynamically_stab+=1./Nsimuls;
            }
          }
          myfile << g << " " << S << " " << dynamically_stab << " ";
        }
      }
      myfile << std::endl;
    }

    myfile.close();

  }catch(error e){
    e.handle();
  }
  return 0;
}
