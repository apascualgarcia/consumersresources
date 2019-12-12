#include "../../include/CRModel.h"
#include <cmath>
ntype distance_between_equilibria(const Extinction& extinct){
  ntype distance_spec = 0., distance_res = 0.;
  size_t NR = extinct.old_Req.size(), NS = extinct.old_Seq.size();
  for(size_t nu=0; nu < NR;++nu){
    ntype to_add = (1.-extinct.new_Req[nu]/extinct.old_Req[nu])/NR;
    distance_res += (to_add*to_add);
  }
  for(size_t i = 0; i < NS; ++i){
    ntype to_add = (1.-extinct.new_Seq[i]/extinct.old_Seq[i])/NS;
    distance_spec+=(to_add*to_add);
  }
  return sqrt(distance_spec+distance_res);
}
