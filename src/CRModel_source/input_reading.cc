#include "../../include/CRModel.h"

gammamode string_to_gamma_mode(std::string mode){
  if(mode=="random_val"){
    return gammamode(random_val);
  }else if(mode=="nested"){
    return gammamode(nested);
  }else if(mode=="antinested"){
    return gammamode(antinested);
  }else{
    std::cerr << "Error, that value of gammamode has not been implemented yet or does not exist"<<std::endl;
    std::cerr << "Continuing with gammamode random_val" << std::endl;
    return gammamode(random_val);
  }
}
taumode string_to_tau_mode(std::string mode){
  if(mode=="tau0"){
    return taumode(tau0);
  }else if(mode=="taualpha"){
    return taumode(taualpha);
  }else{
    std::cerr << "Error, that value of taumode has not been implemented yet or does not exist"<<std::endl;
    std::cerr << "Continuing with taumode tau0" << std::endl;
    return taumode(tau0);
  }
}
alphamode string_to_alpha_mode(std::string mode){
  if(mode=="random_structure"){
    return alphamode(random_structure);
  }else if(mode=="no_release_when_eat"){
    return alphamode(no_release_when_eat);
  }else if(mode=="one_release"){
    return alphamode(one_release);
  }else if(mode=="optimal_matrix"){
    return alphamode(optimal_matrix);
  }else{
    error err("Error, that value of alphamode has not been implemented yet or does not exist.");
    throw err;
  }
}

buildingmode string_to_building_mode(std::string mode){
  if(mode=="use_l"){
    return buildingmode(use_l);
  }else if(mode=="use_m"){
    return buildingmode(use_m);
  }else{
    error err("Error, that value of building mode has not been implemented yet or does not exist.");
    throw err;
  }
}

eqmode string_to_eq_mode(std::string mode){
    if(mode=="one_extinct"){
      return eqmode(oneextinct);
    }else if(mode=="convergence"){
      return eqmode(convergence);
    }else{
      error err("Error, that value of equilibrium_mode has not been implemented yet or does not exist.");
      throw err;
    }
}
