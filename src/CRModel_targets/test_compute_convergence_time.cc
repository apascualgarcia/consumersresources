#include "../../include/CRModel.h"

int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    unsigned int Nsimuls = 1;
    double av_millisec=0.;

    for(size_t i=0; i < Nsimuls; ++i){
      CRModel model(metaparams);
      eqmode equilibre(convergence);
      writemode ecriture(false, std::cout);

      model.perturb_parameters(metaparams.perturb_parameters);
      std::chrono::time_point<std::chrono::system_clock> start, end;
      start = std::chrono::system_clock::now();
      model.evolve_until_equilibrium(metaparams.convergence_threshold, equilibre,ecriture);
      end = std::chrono::system_clock::now();

      int elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
      av_millisec+= double(elapsed_seconds)/Nsimuls;
    }
    std::cout << "Average CPU time for convergence : " << av_millisec << std::endl;



  }catch(error e){
    e.handle();
  }

  return 0;
}
