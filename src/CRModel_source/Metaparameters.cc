#include "../../include/CRModel.h"
#include "../../include/ConfigFile.tpp"

Metaparameters::Metaparameters(int argc, char *argv[]){
  std::string inputPath("configuration.in"); // Fichier d'input par defaut
  if (argc > 1) // Fichier d'input specifie par l'utilisateur ("./Exercice3
                // config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans
                                    // une "map" de strings.

  for (int i(2); i < argc; ++i) // Input complementaires ("./Onde
                                // config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  this->gamma0 = configFile.get<ntype>("gamma_0");
  this->alpha0 = configFile.get<ntype>("alpha_0");
  this->sigma0 = configFile.get<ntype>("sigma_0");
  this->p = configFile.get<ntype>("p");
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
  this->l0 = configFile.get<ntype>("l_0");
  this->budget_constraint = configFile.get<bool>("budget_constraint");
  this->seed_number = configFile.get<unsigned int>("seed_number");
  this->tf = configFile.get<ntype>("tf");
  this->perturb_eq = configFile.get<ntype>("perturbation_equilibrium");
  this->perturb_parameters = configFile.get<ntype>("perturbation_parameters");
  this->equilibrium = string_to_eq_mode(configFile.get<std::string>("equilibrium_mode"));
  if(this->verbose > 0){
    std::cout << "Loading metaparameters from " << inputPath << std::endl;
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

ntype Metaparameters::feasible_alpha_max(ntype alpha_accuracy)const{
  ntype alpha_max =0.;
  gsl_function F;
  Metaparameters meta = *this;

  Solver_Parameters solv;
  solv.metaparameters = &meta;
  solv.Nsimul=1000;
  solv.target=0.95;

  F.function = &function_proba_feasability_solver;
  F.params = &solv;

  /* we first find the alpha_max for which the feasability probability is roughly 0.98 */
  alpha_max = find_zero(&F, this->verbose, interval(0, this->physical_maximum_alpha0()));
  meta.alpha0 = alpha_max;

  /*  we take a tiny step back and wait for the first alpha which is = 1
      this works because we know the shape of the feasability vs alpha curve */
  while(find_feasability_probability(meta) < 1.){
    alpha_max -= alpha_accuracy;
    meta.alpha0=alpha_max;
  }

  if(this->verbose > 0){
    std::cout << "Maximum feasible alpha0 for " << this->foodmatrixpath << " is " << alpha_max << std::endl;
  }

  return alpha_max;
}
