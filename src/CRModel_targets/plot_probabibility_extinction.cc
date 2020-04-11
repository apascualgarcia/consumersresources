#include "../../include/CRModel.h"


int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    nvector delta0_interval = linear_interval(0.02,0.1, 50);
    unsigned int Nsimul=100;
    std::ofstream myfile = open_external_file_truncate(metaparams.save_path);
    myfile << "#" << metaparams.foodmatrixpath << std::endl;
    for(auto delta : delta0_interval){
      ntype proba = probability_of_extinction_greather_than_one(&metaparams, delta, Nsimul);
      myfile << delta << " " << proba << " ";
    }

    myfile << std::endl;
    myfile.close();


  }catch(error e){
    e.handle();
  }
  return 0;
}
