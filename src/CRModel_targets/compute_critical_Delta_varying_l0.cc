#include <iostream>
#include "CRModel.h"
#include <fstream>
using namespace std;

/*  What's different between this file and src/CRModel_targets/compute_critical_Delta_matrices.cc is that
    the present file scans over a set of alpha's for the same matrix, while the other keeps the same
    metaparameters but scans over a set of matrices */

int main(int argc, char * argv[]){
  Metaparameters metaparams(argc, argv);

  vector<ntype> l0_vals;
  ntype lmax = 11092, lmin=400;
  unsigned int Nl0=10;
  for(size_t i=0; i < Nl0; ++i){
    l0_vals.push_back(lmin+(lmax-lmin)/(Nl0-1)*i);
  }

  std::ofstream myfile;
  myfile.open(metaparams.save_path, std::ios::trunc);
  bool save_success(false);
  if(not(myfile.is_open())){
    std::cerr << "Could not open " << metaparams.save_path << " to write the new equilibrium of the system" << std::endl;
  }else{
    if(metaparams.verbose > 0){
      std::cout << "Successfully opened " << metaparams.save_path <<", attempting now to find the critical delta of the matrix for the different l0's " << std::endl;
    }
    save_success = true;
    for(size_t i = 0; i < l0_vals.size();++i){
        metaparams.l0 = l0_vals[i];
        statistics delta = compute_critical_Delta(metaparams);
        std::cout << "Computed critical delta for l0=" << l0_vals[i] << std::endl;
        myfile << l0_vals[i] << " " << delta.mean_ << " " << delta.std_deviation_ << std::endl;
    }
  }
  myfile.close();
  return 0;
}
