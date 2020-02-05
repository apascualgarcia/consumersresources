#include <iostream>
#include <fstream>
#include "CRModel.h"

using namespace std;

/*  This script quantifies how many feasible systems you build for certain parameters. This will allow us
    to increase syntrophy while still building a reasonable amount of systems for the same time */
int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    /* loading the different matrices */
    std::vector<std::string> matrices=load_food_matrix_list(metaparams.foodmatrixpath);
    for(size_t i=0; i < matrices.size();++i){
      std::cout << matrices[i] << std::endl;
    }

    /* we attempt to build a Nsystems systems for each point*/
    unsigned int Nsystems = 100;
    /* we check for N_alpha values of alpha between alpha_min and alpha_max*/
    unsigned int N_alpha = 100;
    ntype alpha_min = 0.;
    ntype alpha_max = 0.2;

    nvector alpha_range;
    for(size_t i=0; i < N_alpha; ++i){
      alpha_range.push_back(alpha_min+(alpha_max-alpha_min)*i/(N_alpha-1));
    }

    std::ofstream myfile;
    myfile.open(metaparams.save_path, std::ofstream::out | std::ofstream::trunc);
    bool save_success(false);
    if(not(myfile.is_open())){
      std::cerr << "Could not open " << metaparams.save_path << " to write the critical feasability probabilities" << std::endl;
    }else{
      if(metaparams.verbose > 0){
        std::cout << "Successfully opened " << metaparams.save_path <<", attempting now to find the critical feasability probability of every listed matrix " << std::endl;
      }
      for(size_t i=0; i < matrices.size();++i){
        metaparams.foodmatrixpath = matrices[i];
        nvector proba_range;
        bool completely_feasible = true;
        myfile << matrices[i] << " ";
        for(size_t j=0; j < N_alpha;++j){
          metaparams.alpha0 = alpha_range[j];
          ntype proba = find_feasability_probability(metaparams, Nsystems);
          if(proba*proba < 1.){
            completely_feasible = false;
          }
          myfile << alpha_range[j] << " " << proba << " ";
          proba_range.push_back(proba);
        }
        myfile << std::endl;
      }
    }
    myfile.close();    
  }catch(error e){
    e.handle();
  }
  return 0;
}
