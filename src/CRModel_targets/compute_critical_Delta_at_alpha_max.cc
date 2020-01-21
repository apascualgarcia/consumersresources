#include <iostream>
#include <fstream>
#include "CRModel.h"

using namespace std;

/*  This script quantifies how many feasible systems you build for certain parameters. This will allow us
    to increase syntrophy while still building a reasonable amount of systems for the same time */
int main(int argc, char * argv[]){
  Metaparameters metaparams(argc, argv);
  unsigned int Npoints = 3;
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

  std::ofstream myfile;
  myfile.open(metaparams.save_path, std::ofstream::out | std::ofstream::trunc);
  bool save_success(false);
  if(not(myfile.is_open())){
    std::cerr << "Could not open " << metaparams.save_path << " to write the critical feasability probabilities" << std::endl;
  }else{
    if(metaparams.verbose > 0){
      std::cout << "Successfully opened " << metaparams.save_path <<", attempting now to find the critical feasability probability of every listed matrix " << std::endl;
    }
    const ntype initial_a0 = metaparams.alpha0;
    for(size_t i=0; i < matrices.size();++i){
      metaparams.foodmatrixpath = matrices[i];
      ntype alpha_max = metaparams.feasible_alpha_max();
      myfile<< matrices[i] << " ";
      for(size_t j=0; j < Npoints; ++j){
        metaparams.alpha0 = alpha_max*j/(Npoints-1);
        statistics critical_delta = compute_critical_Delta(metaparams);
        myfile << metaparams.alpha0 << " " << critical_delta.mean_ << " " << critical_delta.std_deviation_ << " ";
      }
      myfile << std::endl;
    }
  }
  myfile.close();

  return 0;
}
