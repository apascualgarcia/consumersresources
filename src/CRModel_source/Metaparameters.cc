#include "../../include/CRModel.h"
#include "../../include/ConfigFile.tpp"

Metaparameters::Metaparameters(int argc, char *argv[]){
  std::string inputPath("config/configuration.in"); // Fichier d'input par defaut
  if (argc > 1) // Fichier d'input specifie par l'utilisateur ("./Exercice3
                // config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans
                                    // une "map" de strings.

  for (int i(2); i < argc; ++i) // Input complementaires ("./Onde
                                // config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  this->gamma0 = configFile.get<ntype>("gamma0");
  this->alpha0 = configFile.get<ntype>("alpha0");
  this->sigma0 = configFile.get<ntype>("sigma0");
  this->p = configFile.get<ntype>("p");
  this->m0 = configFile.get<ntype>("m0");
  this->R0 = configFile.get<ntype>("R0");
  this->S0 = configFile.get<ntype>("S0");
  this->NR = configFile.get<unsigned int>("NR");
  this->NS = configFile.get<unsigned int>("NS");
  this->gamma_mode = string_to_gamma_mode(configFile.get<std::string>("gamma_mode"));
  this->tau_mode = string_to_tau_mode(configFile.get<std::string>("tau_mode"));
  this->alpha_mode = string_to_alpha_mode(configFile.get<std::string>("alpha_mode"));
  this->foodmatrixpath = configFile.get<std::string>("path_to_food_matrix");
  this->verbose=configFile.get<unsigned int>("verbose-level");
  this->energy_constraint=configFile.get<bool>("energy_constraint");
  this->nb_attempts = configFile.get<unsigned int>("number_of_attempts");
  this->save_path = configFile.get<std::string>("path_to_save_file");
  this->epsilon = configFile.get<ntype>("epsilon");
  this->l0 = configFile.get<ntype>("l0");
  this->budget_constraint = configFile.get<bool>("budget_constraint");
  this->seed_number = configFile.get<unsigned int>("seed_number");
  this->tf = configFile.get<ntype>("tf");
  this->perturb_eq = configFile.get<ntype>("perturbation_equilibrium");
  this->perturb_parameters = configFile.get<ntype>("perturbation_parameters");
  this->equilibrium = string_to_eq_mode(configFile.get<std::string>("equilibrium_mode"));
  this->convergence_threshold = configFile.get<ntype>("convergence_threshold");
  this->building_mode = string_to_building_mode(configFile.get<std::string>("building_mode"));
  if(this->verbose > 0){
    std::cout << "Metaparameters loaded : " << *this << std::endl;
  }
  initialize_random_engine(*this);

  return;
}

ntype Metaparameters::physical_maximum_alpha0() const{
  ntype maxa0=0.;

  maxa0 = (1+this->epsilon)*(1+this->epsilon)/(1-this->epsilon);
  maxa0 *= (1-(1-this->epsilon)*this->sigma0);
  maxa0 *= this->NR*this->gamma0*this->R0;

  return maxa0;
}

ntype Metaparameters::minimum_S0_guaranteed_feasability()const{
  ntype result=(1.-this->epsilon)*this->l0;
  result = result/(this->NS*this->gamma0*this->R0*gsl_pow_3(1+this->epsilon)-this->alpha0*gsl_pow_2(1-this->epsilon));
  return result;
}

ntype Metaparameters::feasible_S0_max(ntype S0_accuracy) const{
  ntype S0_max=0.;
  gsl_function F;
  Metaparameters meta = *this;

  Solver_Parameters solv;
  solv.metaparameters = &meta;
  solv.Nsimul=100;
  solv.target=0.99;

  F.function = &function_proba_feasability_solver_S0;
  F.params = &solv;

  try{
    S0_max = find_zero(&F, this->verbose, interval(0, this->R0));
    meta.S0 = S0_max;

    if(this->verbose > 0){
      std::cout << "Maximum feasible S0 for " << this->foodmatrixpath << " is " << S0_max << std::endl;
    }

    if(is_an_error(S0_max)){
      error err("Could not find an appropriate value for the feasible S0_max. The exception value is returned.",1);
      throw err;
    }
    if(S0_max < 0){
      error err("Found a negative S0_max. The exception value will be returned.", 1);
      throw err;
    }

  }catch(error e){
    e.handle();
    return NUMERICAL_ERROR;
  }

  return S0_max;
}

ntype Metaparameters::feasible_alpha_max(ntype alpha_accuracy)const{
  ntype alpha_max =0.;
  gsl_function F;
  Metaparameters meta = *this;

  Solver_Parameters solv;
  solv.metaparameters = &meta;
  solv.Nsimul=100;
  solv.target=0.99;

  F.function = &function_proba_feasability_solver;
  F.params = &solv;

  /* we first find the alpha_max for which the feasability probability is roughly 0.98 */
  try{
    alpha_max = find_zero(&F, this->verbose, interval(0, this->physical_maximum_alpha0()));
    meta.alpha0 = alpha_max;
    /*  we take a tiny step back and wait for the first alpha which is = 1
        this works because we know the shape of the feasability vs alpha curve */
    /*
    unsigned int steps=0;
    while(find_feasability_probability(meta) < 1. and steps<=1000){
      alpha_max -= alpha_accuracy;
      meta.alpha0=alpha_max;
      steps+=1;
    }
    */

    if(this->verbose > 0){
      std::cout << "Maximum feasible alpha0 for " << this->foodmatrixpath << " is " << alpha_max << std::endl;
    }

    if(is_an_error(alpha_max) || alpha_max < 0){
      error err("Could not find an appropriate value for the feasible alpha_max. The exception value is returned.",1);
      throw err;
    }

  }catch(error e){
    e.handle();
    return NUMERICAL_ERROR;
  }


  return alpha_max;
}

ntype Metaparameters::feasible_alpha_min(ntype alpha_accuracy) const{
  ntype alpha_min=0.;
  ntype upper_bound = this->feasible_alpha_max(alpha_accuracy)-alpha_accuracy;

  gsl_function F;
  Metaparameters meta = *this;

  Solver_Parameters solv;
  solv.metaparameters = &meta;
  solv.Nsimul=1000;
  solv.target=0.95;

  F.function = &function_proba_feasability_solver;
  F.params = &solv;

  /* we first find the alpha_max for which the feasability probability is roughly 0.98 */
  alpha_min = find_zero(&F, this->verbose, interval(0, upper_bound));
  meta.alpha0 = alpha_min;

  if(this->verbose > 0){
    std::cout << "Minimum feasible alpha0 for " << this->foodmatrixpath << " is " << alpha_min << std::endl;
  }



  return alpha_min;
}
