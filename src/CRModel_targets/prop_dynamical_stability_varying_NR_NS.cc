#include "../../include/CRModel.h"
#include <fstream>

int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    unsigned int Nsimuls = 100;
    unsigned int NRmin=1, NRmax=25;
    unsigned int NSmin=1, NSmax=25;

    /* will write a matrix with NR as line index and NS as column index */

    metaparams.foodmatrixpath = "matrices/Nr100_Nc100/RandTrix_Nr100_Nc100_Nest1_Conn1.txt";
    std::ofstream myfile, myfile_eff;
    myfile.open(metaparams.save_path, std::ofstream::out | std::ofstream::trunc);
    myfile_eff.open(metaparams.save_path+"_eff", std::ofstream::out | std::ofstream::trunc);
    bool save_success(false);
    if(not(myfile.is_open())){
      std::cerr << "Could not open " << metaparams.save_path << std::endl;
    }else{
      if(metaparams.verbose > 0){
        std::cout << "Successfully opened " << metaparams.save_path <<", now attempting to compute the phase diagram of the given metaparameters" << std::endl;
      }
      for(size_t i=0; i < NRmax-NRmin+1; ++i){
        unsigned int current_NR = NRmin +(NRmax-NRmin)*i/(NRmax-NRmin);
        for(size_t j=0; j < NSmax-NSmin+1; ++j){
          unsigned int current_NS = NSmin+(NSmax-NSmin)*j/(NSmax-NSmin);
          metaparams.NR = current_NR;
          metaparams.NS = current_NS;
          stability stab_full = compute_proportion_stability(metaparams, Nsimuls, full);
          stability stab_eff = compute_proportion_stability(metaparams, Nsimuls, effective);
          myfile << current_NR << " " << current_NS << " " << stab_full << std::endl;
          myfile_eff << current_NR << " " << current_NS << " " << stab_eff << std::endl;
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
