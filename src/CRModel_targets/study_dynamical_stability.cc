#include "../../include/CRModel.h"

int main(int argc, char * argv[]){
  try{
    Metaparameters metaparams(argc, argv);
    initialize_random_engine(metaparams);
    ntype inter_start = 0., inter_end = 1.;
    size_t Npoints = 10;
    stabilitymode stab_mode = dynamical;
    for(size_t i=0; i < Npoints;++i){
      ntype delta = inter_start+i*(inter_end-inter_start)/(Npoints-1);
      stability_metrics stab_metrics = compute_stability_metrics(metaparams, delta, 100, stab_mode);
    }
  }catch(error e){
    e.handle();
  }

  return 0;
}
