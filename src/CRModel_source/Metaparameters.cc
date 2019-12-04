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
  return;
}
