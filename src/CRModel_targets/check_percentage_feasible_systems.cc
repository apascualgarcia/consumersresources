#include <iostream>
#include <fstream>
#include "CRModel.h"

using namespace std;

/*  This script quantifies how many feasible systems you build for certain parameters. This will allow us
    to increase syntrophy while still building a reasonable amount of systems for the same time */
int main(int argc, char * argv[]){
  Metaparameters metaparams(argc, argv);
  /* loading the different matrices */
  std::vector<std::string> matrices;
  std::ifstream in(metaparams.foodmatrixpath);
  if (!in) {
    std::cerr << "Cannot open file containing the list of matrices " << metaparams.foodmatrixpath << " to measure their critical probability feasability" << std::endl;
  }else{
    do{
      std::string a;
      in >> a;
      matrices.push_back(a);
    }while(!in.eof());
    /* remove last string if white space */
    std::string str = matrices[matrices.size()-1];
    if(str.find_first_not_of(' ') == std::string::npos){
      matrices.pop_back();
    }
  }
  in.close();
  for(size_t i=0; i < matrices.size();++i){
    std::cout << matrices[i] << std::endl;
  }

  /* we attempt to build a Nsystems systems for each point*/
  unsigned int Nsystems = 100;
  /* we check for N_alpha values of alpha between alpha_min and alpha_max*/
  unsigned int N_alpha = 100;
  ntype alpha_min = 4.;
  ntype alpha_max = 16.5;

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
      for(size_t j=0; j < N_alpha;++j){
        metaparams.alpha0 = alpha_range[j];
        ntype proba = find_feasability_probability(metaparams, Nsystems);
        if(proba*proba < 1.){
          completely_feasible = false;
        }
        proba_range.push_back(proba);
      }
      if(not(completely_feasible)){
        myfile << matrices[i] << " ";
        myfile << proba_range << std::endl;
      }else{
        std::cerr << "Matrix " << matrices[i] << " has a critical alpha outside the range considered " << std::endl;
      }
    }
  }
  myfile.close();

  return 0;
}
