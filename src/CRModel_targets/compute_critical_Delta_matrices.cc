#include <iostream>
#include "CRModel.h"
#include <fstream>
using namespace std;

int main(int argc, char * argv[]){
  Metaparameters metaparams(argc, argv);
  initialize_random_engine(metaparams);
  vector<string> matrices_path;
  std::ifstream in(metaparams.foodmatrixpath);
  if (!in) {
    std::cerr << "Cannot open file " << metaparams.foodmatrixpath << " to compute critical deltas" << std::endl;
  }else{
    do{
      std::string a;
      in >> a;
      matrices_path.push_back(a);
    }while(!in.eof());
  }
  in.close();
  std::ofstream myfile;
  myfile.open(metaparams.save_path, std::ios::app);
  if(metaparams.verbose > 0){
    std::cout << "The critical deltas for the following matrices will be computed :";
    for(size_t i = 0 ; i < matrices_path.size(); ++i){
      std::cout << std::endl << matrices_path[i];
    }
  }
  bool save_success(false);
  if(not(myfile.is_open())){
    std::cerr << "Could not open " << metaparams.save_path << " to write the new equilibrium of the system" << std::endl;
  }else{
    if(metaparams.verbose > 0){
      std::cout << "Successfully opened " << metaparams.save_path <<", attempting now to find the critical delta of every listed matrix " << std::endl;
    }
    save_success = true;
    for(size_t i = 0; i < matrices_path.size();++i){
        metaparams.foodmatrixpath = matrices_path[i];
        delta_solver solv_params = {fitmode(sigmoidal),eqmode(oneextinct)};
        statistics delta = compute_critical_Delta(metaparams, solv_params);
        std::cout << "Computed critical delta for " << matrices_path[i] << std::endl;
        myfile << matrices_path[i] << " " << delta.mean_ << " " << delta.std_deviation_ << std::endl;
    }
  }
  myfile.close();
  return 0;
}
