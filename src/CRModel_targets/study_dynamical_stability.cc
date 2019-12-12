#include "../../include/CRModel.h"

int main(int argc, char * argv[]){
  Metaparameters metaparams(argc, argv);
  initialize_random_engine(metaparams);
  ntype inter_start = 0., inter_end = 1.;
  size_t Npoints = 10;
  for(size_t i=0; i < Npoints;++i){
    ntype delta = inter_start+i*(inter_end-inter_start)/(Npoints-1);
    //double proba = probability_of_extinction_greather_than_one(&metaparams, delta, 100, stabilitymode(dynamical));
    ntype dist = average_distance_between_equilibria(&metaparams, delta, 100, stabilitymode(dynamical));
  }

  return 0;
}
