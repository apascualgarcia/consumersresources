#include "CRModel.h"
#include <iostream>
#include <fstream>
#include <gsl/gsl_statistics_double.h>
#include <iomanip>

/*  this script allows to get the average m0 and d0 for
    a given, fixed set of metaparameters */

int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    initialize_random_engine(metaparams);

    /* the number of systems used to estimate m0 and d0 */
    unsigned int Nsimul = 1;

    /* we write the results on the file path given as input */
    std::ofstream myfile;
    myfile.open(metaparams.save_path, std::ios::app);
    if(not(myfile.is_open())){
      std::cerr << "Could not open " << metaparams.save_path << " to write the temporal evolution" << std::endl;
    }else{
      for(size_t i = 0 ; i < Nsimul; ++i){
        myfile <<  std::fixed << std::setprecision(5);
        CRModel model(metaparams);
        double m0 = model.get_m0(), d0 = model.get_d0();
        nvector death_resources = model.get_m();
        double death_res[death_resources.size()];
        for(size_t j = 0; j<death_resources.size();++j){
          death_res[j] = death_resources[j];
        }
        nvector death_consumers = model.get_d();
        double death_cons[death_consumers.size()];
        for(size_t j = 0; j<death_consumers.size();++j){
          death_cons[j] = death_consumers[j];
        }
        std::cout << death_resources << std::endl;
        std::cout << "Mean : " << gsl_stats_mean(death_res, 1, death_resources.size()) << std::endl;
        std::cout << "Std : " <<gsl_stats_sd(death_res,1, death_resources.size()) << std::endl;
        std::cout << death_consumers << std::endl;
        std::cout << "Mean : " << gsl_stats_mean(death_cons, 1, death_consumers.size()) << std::endl;
        std::cout << "Std : " <<gsl_stats_sd(death_cons,1, death_consumers.size()) << std::endl;
        myfile << metaparams.alpha0 << " " << death_resources.size() << " ";
        myfile << death_consumers.size() << " "<<  m0 << " " << d0 << std::endl;
      }
    }
  }catch(error e){
    e.handle();
  }
  return 0;
}
