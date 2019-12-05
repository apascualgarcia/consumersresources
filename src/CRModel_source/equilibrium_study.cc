#include "../../include/CRModel.h"
#include <gsl/gsl_statistics.h>

Extinction_statistics compute_average_extinction(Metaparameters* metaparams, const ntype& Delta, unsigned int Nsimul){
  ntype convergence_threshold = 1e-6;
  Extinction_statistics av_extinct;
  av_extinct.t_eq.mean = 0.;
  av_extinct.t_eq.std_deviation= 0.;
  av_extinct.extinct.mean = 0.;
  av_extinct.extinct.std_deviation = 0.;
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

  av_extinct.t_eq.mean = gsl_stats_mean(teq, 1, Nsimul);
  av_extinct.t_eq.std_deviation = gsl_stats_sd_m(teq, 1, Nsimul, av_extinct.t_eq.mean);

  av_extinct.extinct.mean = gsl_stats_mean(extinctions, 1, Nsimul);
  av_extinct.extinct.std_deviation = gsl_stats_sd_m(extinctions, 1, Nsimul, av_extinct.extinct.mean);

  if(metaparams->verbose > 0){
    std::cout << "Average extinction for Delta = " << Delta << " is " << av_extinct.extinct.mean;
    std::cout << " +/- " << av_extinct.extinct.std_deviation << " (" <<Nsimul << " runs)" << std::endl;
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
  double av_number_extinct = double(ext.extinct.mean);

  return av_number_extinct;
}
double probability_of_extinction_greather_than_one(Metaparameters* metaparams, const ntype& Delta, unsigned int Nsimul){
  ntype convergence_threshold = 1e-6;
  double probability_ext_gtone = 0.;

  if(metaparams->verbose > 0){
    std::cout << "Computing probability of getting one or more extinctions for Delta=" << Delta << std::endl;
  }

  for(size_t i = 0; i < Nsimul; ++i){
    CRModel model(*metaparams);
    model.perturb_parameters(Delta);
    Extinction new_equilib = model.evolve_until_equilibrium(convergence_threshold, eqmode(oneextinct));
    if(new_equilib.extinct >=1){
      probability_ext_gtone+=1./Nsimul;
    }
  }
  if(metaparams->verbose>0){
    std::cout << "Probability of getting one or more extinctions for Delta=" << Delta;
    std::cout << " is " << probability_ext_gtone << std::endl;
  }
  return probability_ext_gtone;
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
