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

  /* we attempt to build Nsystems systems for each point*/
  unsigned int Nsystems = 10000;


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
      myfile << matrices[i] << " ";
      statistics alpha_crit = compute_critical_alpha(metaparams, 1e-4, sigmoidal);
      myfile << alpha_crit.mean_ << " " << alpha_crit.std_deviation_ << std::endl;
    }
  }
  myfile.close();

  return 0;
}
