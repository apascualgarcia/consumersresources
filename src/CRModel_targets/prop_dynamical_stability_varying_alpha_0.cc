#include "../../include/CRModel.h"
#include <fstream>

int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    unsigned int Nsimuls = 100;
    ntype alpha_min=0., alpha_max=metaparams.feasible_alpha_max();
    unsigned int alpha_points=100;

    std::ofstream myfile, myfile_eff;
    myfile.open(metaparams.save_path, std::ofstream::out | std::ofstream::trunc);
    myfile_eff.open(metaparams.save_path+"_eff", std::ofstream::out | std::ofstream::trunc);

    bool save_success(false);
    if(not(myfile.is_open()) or not(myfile_eff.is_open())){
      std::cerr << "Could not open " << metaparams.save_path << std::endl;
    }else{
      if(metaparams.verbose > 0){
        std::cout << "Successfully opened " << metaparams.save_path <<", now attempting to compute the phase diagram of the given metaparameters" << std::endl;
      }
      nvector alpha_interval;
      for(size_t i=0; i < alpha_points; ++i){
        alpha_interval.push_back(alpha_min+(alpha_max-alpha_min)*i/(alpha_points-1));
      }

      for(size_t j=0; j < alpha_points; ++j){
        metaparams.alpha0=alpha_interval[j];
        /* we compute the stability of both full and effective models */
        stability stab_full = compute_proportion_stability(metaparams, Nsimuls, full);
        stability stab_eff = compute_proportion_stability(metaparams, Nsimuls, effective);
        std::cout << "Full stability :" << stab_full << std::endl;
        std::cout << "Effective stability : " << stab_eff << std::endl;
      }

      myfile << std::endl;
    }
    myfile.close();
  }catch(error e){
    e.handle();
  }

  return 0;
}
