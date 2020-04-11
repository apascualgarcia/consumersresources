#include "CRModel.h"

int main(int argc, char* argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    /* we pick a network where we observe patchiness. Take one S0, increase gamma0, measure max eigenvalue */
    nvector gamma0_interval=linear_interval(0.4, 0.5, 1000);
    unsigned int Nsimuls=1000;

    for(auto g0 : gamma0_interval){
      metaparams.gamma0=g0;
      ntype proba_dynamically_stable=0.;
      for(unsigned int i=0; i < Nsimuls; ++i){
        CRModel model(metaparams);
        if(model.is_dynamically_stable()){
          proba_dynamically_stable+=1./Nsimuls;
        }
      }
      std::cout << "Probability of dyn stability for gamma0="<< metaparams.gamma0 << " : " << proba_dynamically_stable << std::endl;
    }


  }catch(error e){
    e.handle();
  }
  return 0;
}
