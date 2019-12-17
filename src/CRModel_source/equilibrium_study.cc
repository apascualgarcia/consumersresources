#include "../../include/CRModel.h"
#include <gsl/gsl_statistics.h>

Extinction_statistics compute_average_extinction(Metaparameters* metaparams, const ntype& Delta, unsigned int Nsimul){
  ntype convergence_threshold = 1e-6;
  Extinction_statistics av_extinct;
  av_extinct.t_eq.mean_ = 0.;
  av_extinct.t_eq.std_deviation_= 0.;
  av_extinct.extinct.mean_ = 0.;
  av_extinct.extinct.std_deviation_ = 0.;
  av_extinct.new_Req.means = nvector(metaparams->NR, 0.);
  av_extinct.new_Seq.means = nvector(metaparams->NS, 0.);

  foodmatrix food_matrix = load_food_matrix(*metaparams);
  if(metaparams->verbose > 0){
    std::cout << "Computing average extinction for delta="<<Delta << " with food matrix " << metaparams->foodmatrixpath << std::endl;
  }

  double teq[Nsimul];
  double extinctions[Nsimul];

  for(size_t i = 0; i < Nsimul; ++i){
    CRModel model(food_matrix,*metaparams);
    model.perturb_parameters(Delta);
    Extinction new_equilib = model.evolve_until_equilibrium(convergence_threshold);

    teq[i] = new_equilib.t_eq;
    extinctions[i] = new_equilib.extinct;

    for(size_t j = 0 ; j < metaparams->NS; ++j){
      av_extinct.new_Seq.means[j] += (new_equilib.new_Seq[j]/Nsimul);
    }
    for(size_t mu = 0; mu < metaparams->NR; ++mu){
      av_extinct.new_Req.means[mu] += (new_equilib.new_Req[mu]/Nsimul);
    }
  }

  av_extinct.t_eq.mean_ = gsl_stats_mean(teq, 1, Nsimul);
  av_extinct.t_eq.std_deviation_ = gsl_stats_sd_m(teq, 1, Nsimul, av_extinct.t_eq.mean_);

  av_extinct.extinct.mean_ = gsl_stats_mean(extinctions, 1, Nsimul);
  av_extinct.extinct.std_deviation_ = gsl_stats_sd_m(extinctions, 1, Nsimul, av_extinct.extinct.mean_);

  if(metaparams->verbose > 0){
    std::cout << "Average extinction for Delta = " << Delta << " is ";
    std::cout << av_extinct.extinct<< " (" <<Nsimul << " runs)" << std::endl;
  }
  return av_extinct;
}
double average_number_of_extinctions(double delta, Metaparameters* m, unsigned int Nsimul){
  /* old version
  Solver_Parameters* s = (Solver_Parameters*) params;
  Metaparameters* m = s->metaparameters;
  unsigned int Nsimul = s->Nsimul;
  */

  Extinction_statistics ext = compute_average_extinction(m, ntype(delta), Nsimul);
  double av_number_extinct = double(ext.extinct.mean_);

  return av_number_extinct;
}
double probability_of_extinction_greather_than_one(Metaparameters* metaparams, const ntype& Delta, unsigned int Nsimul, stabilitymode stab_mode){
  ntype convergence_threshold = 1e-9;
  double probability_ext_gtone = 0.;
  std::string stability;

  if(metaparams->verbose > 0){
    std::cout << "Computing probability of getting one or more extinctions for Delta=" << Delta << std::endl;
  }

  switch(stab_mode){
    case structural:{
      stability = "structural";
      for(size_t i = 0; i < Nsimul; ++i){
        CRModel model(*metaparams);
        model.perturb_parameters(Delta);
        Extinction new_equilib = model.evolve_until_equilibrium(convergence_threshold, eqmode(oneextinct));
        if(new_equilib.extinct >=1){
          probability_ext_gtone+=1./Nsimul;
        }
      }
      break;
    };
    case dynamical:{
      stability="dynamical";
      for(size_t i=0; i < Nsimul;++i){
        CRModel model(*metaparams);
        Extinction new_equilib = model.evolve_until_equilibrium_from_abundances(model.perturb_abundances(Delta), convergence_threshold);
        if(new_equilib.extinct >=1){
          probability_ext_gtone += 1./Nsimul;
        }
      }
      break;
    };
    default:{
      std::cerr << "Inexistent type of stability in probability_of_extinction_greather_than_one" << std::endl;
      std::cerr << "Aborting simulation now "<< std::endl;
      abort();
      break;
    }
  }

  if(metaparams->verbose>0){
    std::cout << "Probability of getting one or more extinctions for Delta=" << Delta;
    std::cout << " is " << probability_ext_gtone << std::endl;
  }
  return probability_ext_gtone;
}
statistics distance_between_equilibria(Metaparameters* metaparams, const ntype& delta, unsigned int Nsimul, stabilitymode stab_mode){
    statistics av_dist;
    nvector dist;
    if(metaparams->verbose>0){
      std::cout << "Computing now the average distance between equilibria for delta=" << delta;
    }
    switch(stab_mode){
      case dynamical:{
        if(metaparams->verbose>0){
          std::cout << " for dynamical stability" << std::endl;
        }
        for(size_t i=0; i < Nsimul;++i){
          CRModel model(*metaparams);
          Extinction new_equilib = model.evolve_until_equilibrium_from_abundances(model.perturb_abundances(delta));
          dist.push_back(distance_between_equilibria(new_equilib)/Nsimul);
        }
        av_dist.mean_ = mean(dist);
        av_dist.std_deviation_ = standard_dev(dist);
        break;
      }
      default:{
        std::cerr << "average_distance_between_equilibria not yet implemented for this type of stability mode " << std::endl;
        std::cerr << "Aborting simulation now" << std::endl;
        abort();
        break;
      }
    }
    if(metaparams->verbose>0){
      std::cout << "Average distance between equilibria is " << av_dist << std::endl;
    }

    return av_dist;
}

statistics angle_between_equilibria(Metaparameters* metaparams, const ntype& delta, unsigned int Nsimul, stabilitymode stab_mode){
  statistics angle;
  nvector angles;
  if(metaparams->verbose>0){
    std::cout << "Computing now the average angle between equilibria for delta=" << delta;
  }
  switch(stab_mode){
    case dynamical:{
      if(metaparams->verbose>0){
        std::cout << " for dynamical stability" << std::endl;
      }
      for(size_t i=0; i < Nsimul;++i){
        CRModel model(*metaparams);
        Extinction new_equilib = model.evolve_until_equilibrium_from_abundances(model.perturb_abundances(delta));
        angles.push_back(angle_between_equilibria(new_equilib)/Nsimul);
      }
      angle.mean_ = mean(angles);
      angle.std_deviation_ = standard_dev(angles);
      break;
    }
    default:{
      std::cerr << "average_distance_between_equilibria not yet implemented for this type of stability mode " << std::endl;
      std::cerr << "Aborting simulation now" << std::endl;
      abort();
      break;
    }
  }
  if(metaparams->verbose>0){
    std::cout << "Average angle between equilibria is " << angle << std::endl;
  }

  return angle;
}

double can_find_one_extinction(Metaparameters* metaparams, const ntype& Delta, unsigned int Nsimul){
  ntype convergence_threshold = 1e-6;
  if(metaparams->verbose > 0){
    std::cout << "We check if within " << Nsimul << " runs, we can find one of them where an extinction occurs ";
    std::cout << "after a perturbation of Delta=" << Delta << std::endl;
  }
  for(size_t i=0; i < Nsimul; ++i){
    CRModel model(*metaparams);
    model.perturb_parameters(Delta);
    Extinction new_equilib = model.evolve_until_equilibrium(convergence_threshold, eqmode(oneextinct));
    if(new_equilib.extinct >= 1){
      return 0.;
    }
  }
  return 1.;
}
double can_find_zero_extinction(Metaparameters* metaparams, const ntype& Delta, unsigned int Nsimul){
  ntype convergence_threshold = 1e-6;
  if(metaparams->verbose > 0){
    std::cout << "We check if within " << Nsimul << " runs, we can find one of them where 0 extinction occurs ";
    std::cout << "after a perturbation of Delta=" << Delta << std::endl;
  }
  for(size_t i = 0; i < Nsimul; ++i){
    CRModel model(*metaparams);
    model.perturb_parameters(Delta);
    Extinction new_equilib = model.evolve_until_equilibrium(convergence_threshold, eqmode(oneextinct));
    if(new_equilib.extinct==0){
      return 0.;
    }
  }
  return 1.;
}