#include <iostream>
#include "CRModel.h"

using namespace std;

int main(int argc, char * argv[]){
  Metaparameters metaparams(argc, argv);
  initialize_random_engine(metaparams);
  compute_critical_Delta(metaparams,0.001);
/*
  Extinction av_extinct = compute_average_extinction(&metaparams,metaparams.perturb_parameters,100);
  std::cout << "For a Delta of " << metaparams.perturb_parameters << ", we find an average time to reach a new equilibrium of ";
  std::cout << av_extinct.t_eq << " and on average " << av_extinct.extinct << " species have perished" << std::endl;
  std::cout << "The final average distributions for the equilibria are : " << std::endl;
  std::cout << av_extinct.new_Req << std::endl << av_extinct.new_Seq << std::endl;
*/
  return 0;
}
