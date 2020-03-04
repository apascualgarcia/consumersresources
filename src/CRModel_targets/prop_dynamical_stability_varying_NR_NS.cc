#include "../../include/CRModel.h"
#include <fstream>

int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    unsigned int Nsimuls = 2000;
    unsigned int NRmin=1, NRmax=25;
    unsigned int NSmin=1, NSmax=25;

    /* will write a matrix with NR as line index and NS as column index */
    metaparams.foodmatrixpath="matrices/Nr100_Nc100/RandTrix_Nr100_Nc100_Nest1_Conn1.txt";
    std::ofstream myfile=open_external_file_truncate(metaparams.save_path);
    std::ofstream myfile_eff=open_external_file_truncate(metaparams.save_path+"_eff");

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
    myfile.close();
    myfile_eff.close();
  }catch(error e){
    e.handle();
  }
  return 0;
}
