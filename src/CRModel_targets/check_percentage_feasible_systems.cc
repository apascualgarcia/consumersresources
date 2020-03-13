#include <iostream>
#include <fstream>
#include "CRModel.h"

using namespace std;

/*  This script quantifies how many feasible systems you build for certain parameters. */
int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    std::vector<std::string> matrix_list = load_food_matrix_list(metaparams.foodmatrixpath);

    ntype min_gamma0=0.01, max_gamma0=1.;
    ntype min_S0=0.01, max_S0=1.;

    nvector gamma0_interval=linear_interval(min_gamma0, max_gamma0, 100);
    nvector S0_interval=linear_interval(min_S0, max_S0, 100);

    std::ofstream myfile = open_external_file_append(metaparams.save_path);
    for(auto mat: matrix_list){
      metaparams.foodmatrixpath=mat;
      myfile << metaparams.foodmatrixpath << " ";
      for(auto g : gamma0_interval){
        for(auto S : S0_interval){
          metaparams.gamma0=g;
          metaparams.S0=S;
          ntype feasability_proba = find_feasability_probability(metaparams);
          myfile << g << " " << S << " " << feasability_proba << " ";
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
