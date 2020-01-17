#include "CRModel.h"
#include <iostream>
#include <fstream>
#include <gsl/gsl_statistics.h>
#include <iomanip>

/*  this script allows to get the average m0 and d0 for
    a given, fixed set of metaparameters */

int main(int argc, char* argv[]){
  Metaparameters metaparams(argc, argv);
  initialize_random_engine(metaparams);

  /* the number of systems used to estimate m0 and d0 */
  unsigned int Nsimul = 1000;

  /* we write the results on the file path given as input */
  std::ofstream myfile;
  myfile.open(metaparams.save_path, std::ios::app);
  if(not(myfile.is_open())){
    std::cerr << "Could not open " << metaparams.save_path << " to write the temporal evolution" << std::endl;
  }else{
    for(size_t i = 0 ; i < Nsimul; ++i){
      myfile <<  std::fixed << std::setprecision(5);
      CRModel model(metaparams);
      nvector death_resources = model.get_m();
      nvector death_consumers = model.get_d();
      myfile << metaparams.alpha0 << " " << death_resources.size() << " ";
      myfile << death_consumers.size() << " "<<  death_resources << " " << death_consumers << std::endl;
    }
  }
  return 0;
}
