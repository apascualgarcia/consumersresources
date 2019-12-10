#include <iostream>
#include "CRModel.h"
#include <fstream>
using namespace std;

/*  What's different between this file and src/CRModel_targets/compute_critical_Delta_matrices.cc is that
    the present file scans over a set of alpha's for the same matrix, while the other keeps the same
    metaparameters but scans over a set of matrices */

int main(int argc, char * argv[]){
  Metaparameters metaparams(argc, argv);
  initialize_random_engine(metaparams);

  vector<ntype> alpha_vals;
  std::ifstream in(metaparams.foodmatrixpath);
  if (!in) {
    std::cerr << "Cannot open file containing the alpha's " << metaparams.foodmatrixpath << " to compute critical deltas" << std::endl;
  }else{
    do{
      ntype a;
      in >> a;
      alpha_vals.push_back(a);
    }while(!in.eof());
  }
  in.close();
  std::ofstream myfile;
  myfile.open(metaparams.save_path, std::ios::app);
  bool save_success(false);
  if(not(myfile.is_open())){
    std::cerr << "Could not open " << metaparams.save_path << " to write the new equilibrium of the system" << std::endl;
  }else{
    if(metaparams.verbose > 0){
      std::cout << "Successfully opened " << metaparams.save_path <<", attempting now to find the critical delta of every listed matrix " << std::endl;
    }
    save_success = true;
    for(size_t i = 0; i < alpha_vals.size();++i){
        metaparams.alpha0 = alpha_vals[i];
        statistics delta = compute_critical_Delta(metaparams, 0.);
        std::cout << "Computed critical delta for alpha0=" << alpha_vals[i] << std::endl;
        myfile << alpha_vals[i] << " " << delta.mean << " " << delta.std_deviation << std::endl;
    }
  }
  myfile.close();
  return 0;
}
