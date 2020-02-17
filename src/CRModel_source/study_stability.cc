#include "../../include/CRModel.h"
#include <gsl/gsl_statistics.h>

stability_metrics compute_stability_metrics(Metaparameters& m, const ntype& delta, unsigned int Nsimul, stabilitymode stab_mode){
  nvector resiliences;
  nvector angles;
  nvector distances;
  nvector extinctions;
  if(m.verbose > 0){
    std::cout << "Now computing the different "<< stab_mode << " stability metrics for delta=" << delta << " ("<<Nsimul<<" runs)" <<std::endl;
  }

  for(size_t i = 0; i < Nsimul; ++i){
    CRModel model(m);
    Extinction new_equilib;
    switch(stab_mode){
      case dynamical:{
        new_equilib = model.evolve_until_equilibrium_from_abundances(model.perturb_abundances(delta));
        break;
      }
      case structural:{
        model.perturb_parameters(delta);
        new_equilib = model.evolve_until_equilibrium(1e-9, convergence);
        break;
      }
    }
    resiliences.push_back(new_equilib.t_eq);
    angles.push_back(angle_between_equilibria(new_equilib));
    distances.push_back(distance_between_equilibria(new_equilib));
    extinctions.push_back(new_equilib.extinct);
    if(m.verbose > 1){
      std::cout << "\t For this run, we found : ";
      std::cout << "resilience = " << resiliences[i] << "; angle = " << angles[i] << "; distance = " << distances[i];
      std::cout << "; number of extinctions = " << extinctions[i] << std::endl;
    }
  }

  statistics resilience(resiliences);
  statistics angle(angles);
  statistics distance(distances);
  statistics extinction(extinctions);

  stability_metrics stab_metr;
  stab_metr.resilience=resilience;
  stab_metr.angle_between_equilibria=angle;
  stab_metr.distance_between_equilibria=distance;
  stab_metr.extinctions=extinction;


  if(m.verbose > 0){
    std::cout << "Overall, we found the following stability metrics : " << stab_metr << std::endl;
  }


  return stab_metr;
}

stability compute_proportion_stability(Metaparameters& meta, unsigned int Nsimuls, CRModelType model_type){
  stability stab_prop = {0., 0., 0.};
  ntype toadd = 1./Nsimuls;
  switch(model_type){
    case full:{
      for(size_t i=0; i < Nsimuls;++i){
        CRModel model(meta);
        switch(model.assess_dynamical_stability()){
          case stable :{
            stab_prop.stable += toadd;
            break;
          }
          case marginal:{
            stab_prop.marginally_stable += toadd;
            break;
          }
          case unstable:{
            stab_prop.unstable += toadd;
            break;
          }
        }
      }
      break;
    }
    case effective:{
      for(size_t i=0; i < Nsimuls;++i){
        EffectiveCRModel model(meta);
        switch(model.assess_dynamical_stability()){
          case stable :{
            stab_prop.stable += toadd;
            break;
          }
          case marginal:{
            stab_prop.marginally_stable += toadd;
            break;
          }
          case unstable:{
            stab_prop.unstable += toadd;
            break;
          }
        }
      }
      break;
    }
  }
  return stab_prop;
}
